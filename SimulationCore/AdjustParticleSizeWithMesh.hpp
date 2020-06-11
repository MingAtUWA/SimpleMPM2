#ifndef __Adjust_Particle_Size_With_Mesh_hpp__
#define __Adjust_Particle_Size_With_Mesh_hpp__

#include "ItemBuffer.hpp"
#include "ItemStack.hpp"
#include "TriangleMeshToParticles.h"

// Important assumptions:
// 
template <typename TriMesh>
class AdjustParticleSizeWithMesh
{
protected:
	typedef TriangleMeshToParticles::Particle Particle;
	typedef typename TriMesh::Element Element;
	typedef typename TriMesh::Node Node;

	TriangleMeshToParticles& mh2ps;

	struct ElemInfo
	{
		size_t id;
		double pcl_vol;
		size_t pcl_num;
		Particle* pcls;
		inline void init()
		{
			pcl_vol = 0.0;
			pcl_num = 0;
			pcls = nullptr;
		}
		inline void add_pcl(Particle* pcl)
		{
			pcl_vol += pcl->vol;
			++pcl_num;
			pcl->next2 = pcls;
			pcls = pcl;
		}
	};
	size_t elem_num;
	ElemInfo *elems_info;

	inline void clear_elems_info()
	{
		if (elems_info)
		{
			delete[] elems_info;
			elems_info = nullptr;
		}
		elem_num = 0;
	}

	inline void alloc_elem_info(size_t num)
	{
		clear_elems_info();
		if (!num) return;
		elem_num = num;
		elems_info = new ElemInfo[elem_num];
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			ElemInfo& ei = elems_info[e_id];
			ei.id = e_id;
			ei.init();
		}
	}

	TriMesh* pmh;
	Element* elems;
	size_t node_num;
	Node* nodes;

public:
	AdjustParticleSizeWithMesh(TriangleMeshToParticles &_mh2ps) :
		mh2ps(_mh2ps), elem_num(0), elems_info(nullptr) {}
	~AdjustParticleSizeWithMesh()
	{
		clear_elems_info();
	}

	void adjust_particles(
		TriMesh& mh,
		size_t min_pcl_num_in_elem,
		double hx, double hy
		)
	{
		alloc_elem_info(mh.get_elem_num());
		elems = mh.get_elems();
		node_num = mh.get_node_num();
		nodes = mh.get_nodes();
		pmh = &mh;

		locate_each_pcl_in_which_elem(mh);

		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			ElemInfo& ei = elems_info[e_id];
			Element& e = elems[e_id];

			if (ei.pcl_num == 0) continue;
			
			if (ei.pcl_num < min_pcl_num_in_elem)
			{
				regenerate_pcls(ei, e, hx, hy, min_pcl_num_in_elem);
			}
			
			adjust_pcl_size(ei, e);
		}

		// debug
		//debug_output_elem_info("test_adj_pcl_size.csv");
	}

	void adjust_particles2(
		TriMesh& mh
		)
	{
		alloc_elem_info(mh.get_elem_num());
		elems = mh.get_elems();
		node_num = mh.get_node_num();
		nodes = mh.get_nodes();
		pmh = &mh;

		locate_each_pcl_in_which_elem(mh);

		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			ElemInfo& ei = elems_info[e_id];
			Element& e = elems[e_id];

			if (ei.pcl_num == 0) continue;

			adjust_pcl_size(ei, e);
		}

		// debug
		//debug_output_elem_info("test_adj_pcl_size.csv");
	}

protected:
	// locate each particle in which element
	// del pcl that is not in mesh
	void locate_each_pcl_in_which_elem(TriMesh& mh)
	{
		Particle *pcl, *pcl2;
		Element* pe;
		pcl = mh2ps.first();
		while (mh2ps.not_end_yet(pcl))
		{
			pe = mh.find_in_which_element(pcl->x, pcl->y);
			if (!pe)
			{
				pcl2 = pcl;
				pcl = mh2ps.next(pcl);
				mh2ps.del_pcl(*pcl2);
				continue;
			}
			elems_info[pe->id].add_pcl(pcl);
			pcl = mh2ps.next(pcl);
		}
	}

	void adjust_pcl_size(ElemInfo &ei, Element &e)
	{
		double vol_ratio = e.area / ei.pcl_vol;
		for (Particle *pcl = ei.pcls; pcl; pcl = pcl->next2)
			pcl->vol *= vol_ratio;
	}

	// for debug
	void debug_output_elem_info(const char* filename)
	{
		std::fstream test_file;
		test_file.open(filename, std::ios::out | std::ios::binary);
		test_file << "elem_id, elem_area, pcl_old_vol, old_error, pcl_new_vol, new_error,\n";
		double pcl_vol;
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			ElemInfo& ei = elems_info[e_id];
			if (ei.pcls)
			{
				pcl_vol = 0.0;
				for (Particle *pcl = ei.pcls; pcl; pcl = pcl->next2)
					pcl_vol += pcl->vol;
				double ea = elems[e_id].area;
				test_file << e_id << ", " << ea << ", "
					<< ei.pcl_vol << ", " << 100.0 * (ei.pcl_vol - ea) / ea << "%,"
					<< pcl_vol << ", " << 100.0 * (pcl_vol - ea) / ea << "%,\n";
			}
		}
		test_file.close();
	}

	// regenerate particles if too little particles in element 
	void regenerate_pcls(
		ElemInfo &ei, Element &e,
		double hx_ori, double hy_ori,
		size_t min_pcl_num_in_elem)
	{
		// get element bounding box
		Rect bbox;
		get_elem_bbox(e, bbox);
		double xlen = bbox.xu - bbox.xl;
		double ylen = bbox.yu - bbox.yl;
		size_t x_num = size_t(ceil(xlen / hx_ori));
		size_t y_num = size_t(ceil(ylen / hy_ori));
		double hx = xlen / double(x_num);
		double hy = ylen / double(y_num);

		// del pcl
		while (true)
		{
			del_pcls_in_elem(ei);
			gen_pcls_in_elem(ei, e,
				bbox.xl, bbox.yl,
				hx, hy, 
				x_num, y_num);
			if (ei.pcl_num > min_pcl_num_in_elem)
				break;
			reselect_x_y_num_hx_and_hy(x_num, y_num, hx, hy, xlen, ylen);
		}
	}

	void get_elem_bbox(Element &e, Rect &bbox)
	{
		Node &n1 = nodes[e.n1];
		Node &n2 = nodes[e.n2];
		Node &n3 = nodes[e.n3];
		bbox.xl = n1.x;
		if (bbox.xl > n2.x)
			bbox.xl = n2.x;
		if (bbox.xl > n3.x)
			bbox.xl = n3.x;
		bbox.xu = n1.x;
		if (bbox.xu < n2.x)
			bbox.xu = n2.x;
		if (bbox.xu < n3.x)
			bbox.xu = n3.x;
		bbox.yl = n1.y;
		if (bbox.yl > n2.y)
			bbox.yl = n2.y;
		if (bbox.yl > n3.y)
			bbox.yl = n3.y;
		bbox.yu = n1.y;
		if (bbox.yu < n2.y)
			bbox.yu = n2.y;
		if (bbox.yu < n3.y)
			bbox.yu = n3.y;
	}

	void del_pcls_in_elem(ElemInfo &ei)
	{
		Particle *pcl;
		while (ei.pcls)
		{
			pcl = ei.pcls;
			ei.pcls = ei.pcls->next2;
			mh2ps.del_pcl(*pcl);
		}
		ei.init();
	}

	void reselect_x_y_num_hx_and_hy(
		size_t &x_num, size_t &y_num,
		double &hx, double &hy,
		double xlen, double ylen)
	{
		// Method 1: increase hx or hy
		double hx_new = xlen / double(x_num + 1);
		double hy_new = ylen / double(y_num + 1);
		if (((hx - hx_new) / hx) > ((hy - hy_new) / hy))
		{
			++x_num;
			hx = hx_new;
		}
		else
		{
			++y_num;
			hy = hy_new;
		}

		// Method 2: increase hx and hy simultaneously
		//++x_num;
		//hx = xlen / double(x_num);
		//++y_num;
		//hy = ylen / double(y_num);
	}
	
	void gen_pcls_in_elem(
		ElemInfo &ei, Element &e,
		double xl, double yl,
		double hx, double hy,
		size_t x_num, size_t y_num
		)
	{
		Particle pcl;
		Particle* ppcl;
		pcl.vol = hx * hy;
		for (size_t y_id = 0; y_id < y_num; ++y_id)
			for (size_t x_id = 0; x_id < x_num; ++x_id)
			{
				pcl.x = xl + (double(x_id) + 0.5) * hx;
				pcl.y = yl + (double(y_id) + 0.5) * hy;
				if (pmh->is_in_triangle(e, pcl.x, pcl.y))
				{
					ppcl = mh2ps.add_pcl(pcl);
					ei.add_pcl(ppcl);
				}
			}


	}
};

#endif