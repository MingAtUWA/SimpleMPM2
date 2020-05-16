#ifndef __Adjust_Particles_With_Triangle_Mesh_hpp__
#define __Adjust_Particles_With_Triangle_Mesh_hpp__

#include "ItemBuffer.hpp"
#include "ItemStack.hpp"
#include "TriangleMeshToParticles.h"

template <typename TriMesh>
class AdjustParticlesWithTriangleMesh
{
protected:
	typedef TriangleMeshToParticles::Particle Particle;

	// for splitting particles
	struct PclPart
	{
		double xl, yl;
		double xu, yu;
		double xc, yc;
		typename TriMesh::Element *pe;
		int depth;
		PclPart *next;
		PclPart &operator= (PclPart &other)
		{
			xl = other.xl;
			yl = other.yl;
			xu = other.xu;
			yu = other.yu;
			xc = other.xc;
			yc = other.yc;
			pe = other.pe;
			depth = other.depth;
			return *this;
		}
	};

	struct ElemPclInfo
	{
		size_t id;
		// for splitting particles
		double pcl_x;
		double pcl_y;
		double pcl_vol;
		bool is_used;
		ElemPclInfo *next;
		// for merging particles
		double area; // merge cloest particles
		Particle *ori_pcls;
		size_t ori_pcl_num;
		Particle *new_pcls;
		Particle *spt_pcls;
		inline void add_ori_pcl(Particle &pcl)
		{
			pcl.next2 = ori_pcls;
			ori_pcls = &pcl;
			++ori_pcl_num;
		}
		inline void add_new_pcl(Particle &pcl)
		{
			pcl.next2 = new_pcls;
			new_pcls = &pcl;
		}
		inline void add_spt_pcl(Particle &pcl)
		{
			pcl.next2 = spt_pcls;
			spt_pcls = &pcl;
		}
		inline void reset()
		{
			ori_pcls = nullptr;
			ori_pcl_num = 0;
			new_pcls = nullptr;
			spt_pcls = nullptr;
		}
	};
	size_t elem_num;
	ElemPclInfo *elem_pcl_infos;
	ElemPclInfo *elem_pcl_info_top; // used elem pcl info
	size_t cur_elem_id; // element where current particle locates
	
	MemoryUtilities::ItemStack<PclPart> pp_stack;
	MemoryUtilities::ItemBuffer<Particle> pcl_buf;

	// helper functions
	void init_elem_pcl_info(TriMesh &mh)
	{
		clear_elem_pcl_info();
		elem_pcl_info_top = nullptr;
		elem_num = mh.get_elem_num();
		typename TriMesh::Element *elems = mh.get_elems();
		elem_pcl_infos = new ElemPclInfo[elem_num];
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			ElemPclInfo &epi = elem_pcl_infos[e_id];
			epi.id = elems[e_id].id;
			// for particle splitting
			epi.is_used = false;
			epi.pcl_vol = 0.0;
			epi.pcl_x = 0.0;
			epi.pcl_y = 0.0;
			// for particle merging
			epi.area = elems[e_id].area;
			epi.reset();
		}
	}
	void clear_elem_pcl_info()
	{
		if (elem_pcl_infos)
		{
			delete[] elem_pcl_infos;
			elem_pcl_infos = nullptr;
		}
		elem_num = 0;
	}
	inline void push_elem_pcl_info(ElemPclInfo &epi)
	{
		if (!epi.is_used)
		{
			epi.is_used = true;
			epi.next = elem_pcl_info_top;
			elem_pcl_info_top = &epi;
		}
	}
	inline ElemPclInfo *pop_elem_pcl_info()
	{
		ElemPclInfo *res = elem_pcl_info_top;
		if (res)
		{
			elem_pcl_info_top = elem_pcl_info_top->next;
			res->is_used = false;
			res->next = nullptr;
		}
		return res;
	}

	inline void add_pcl_prt_to_elem(size_t id, Particle &pcl)
	{
		ElemPclInfo &epi = elem_pcl_infos[id];
		push_elem_pcl_info(epi);

		epi.pcl_x += pcl.x * pcl.vol;
		epi.pcl_y += pcl.y * pcl.vol;
		epi.pcl_vol += pcl.vol;
	}

	inline void merge_pcl_prt(double min_pcl_area)
	{
		ElemPclInfo *pepi;
		while ((pepi = pop_elem_pcl_info()) != nullptr)
		{
			Particle &pcl = *pcl_buf.alloc();
			pcl.x = pepi->pcl_x / pepi->pcl_vol;
			pcl.y = pepi->pcl_y / pepi->pcl_vol;
			pcl.vol = pepi->pcl_vol;
			if (pepi->id == cur_elem_id || pcl.vol > min_pcl_area)
				pepi->add_new_pcl(pcl);
			else
				pepi->add_spt_pcl(pcl);
			// reset
			pepi->pcl_x = 0.0;
			pepi->pcl_y = 0.0;
			pepi->pcl_vol = 0.0;
		}
	}

	void merge_particle(double min_dist_ratio)
	{
		Particle pcl;
		Particle *ppcl1, *ppcl2;
		Particle *prev_pcl, *cur_pcl;
		double min_dist, dist;
		double spt_pcl_vol;

		min_dist_ratio = min_dist_ratio * min_dist_ratio;
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			ElemPclInfo &epi = elem_pcl_infos[e_id];

			// merge pcls close to each other
			min_dist = epi.area * min_dist_ratio;
			while (epi.new_pcls)
			{
				ppcl1 = epi.new_pcls;
				epi.new_pcls = epi.new_pcls->next2;

				pcl.x = ppcl1->x * ppcl1->vol;
				pcl.y = ppcl1->y * ppcl1->vol;
				pcl.vol = ppcl1->vol;

				prev_pcl = nullptr;
				cur_pcl = epi.new_pcls;
				while (cur_pcl)
				{
					dist = (ppcl1->x - cur_pcl->x) * (ppcl1->x - cur_pcl->x)
						 + (ppcl1->y - cur_pcl->y) * (ppcl1->y - cur_pcl->y);
					if (dist < min_dist)
					{
						pcl.x += cur_pcl->x * cur_pcl->vol;
						pcl.y += cur_pcl->y * cur_pcl->vol;
						pcl.vol += cur_pcl->vol;
						// delete this pcl
						if (prev_pcl)
							prev_pcl->next2 = cur_pcl->next2;
						else
							epi.new_pcls = epi.new_pcls->next2;

						ppcl2 = cur_pcl;
						cur_pcl = cur_pcl->next2;
						pcl_buf.del(ppcl2);
						continue;
					}

					prev_pcl = cur_pcl;
					cur_pcl = cur_pcl->next2;
				}
				pcl_buf.del(ppcl1);

				pcl.x /= pcl.vol;
				pcl.y /= pcl.vol;
				ppcl2 = mh2ps.add_pcl(pcl);
				epi.add_ori_pcl(*ppcl2);
			}
			
			// merge spt particle into ori pcls
			if (!epi.spt_pcls)
				continue;
			if (epi.ori_pcls)
			{
				// allocate area of spt particles to ori particles evenly
				spt_pcl_vol = 0.0;
				for (Particle *pcl = epi.spt_pcls; pcl; pcl = pcl->next2)
					spt_pcl_vol += pcl->vol;
				spt_pcl_vol /= double(epi.ori_pcl_num);
				for (Particle *pcl = epi.ori_pcls; pcl; pcl = pcl->next2)
					pcl->vol += spt_pcl_vol;
			}
			else
			{
				// lump into one, shouldn't appear very often
				pcl.x = 0.0;
				pcl.y = 0.0;
				pcl.vol = 0.0;
				for (ppcl1 = epi.spt_pcls; ppcl1; ppcl1 = ppcl1->next2)
				{
					pcl.x += ppcl1->x * ppcl1->vol;
					pcl.y += ppcl1->y * ppcl1->vol;
					pcl.vol += ppcl1->vol;
				}
				pcl.x /= pcl.vol;
				pcl.y /= pcl.vol;
				mh2ps.add_pcl(pcl);
			}
		}
	}
	
	TriangleMeshToParticles &mh2ps;
	
public: 
	AdjustParticlesWithTriangleMesh(TriangleMeshToParticles &_mh2ps) :
		mh2ps(_mh2ps), elem_pcl_infos(nullptr)
	{
		pp_stack.set_page_size(200);
		pcl_buf.set_page_size(200);
	}
	~AdjustParticlesWithTriangleMesh() {}
	
	// particle whose area < min_pcl_area, will merge it to other particles
	void distribute_points_area_to_mesh(
		TriMesh &mh, double pcl_hx, double pcl_hy,
		double min_area_ratio = 0.0,
		double min_len_ratio = 0.0,
		int max_depth = 5)
	{
		if (max_depth <= 0) return;

		++max_depth;
		double *pcl_vols = new double[max_depth];
		pcl_vols[0] = pcl_hx * pcl_hy;
		for (size_t d_id = 1; d_id < max_depth; ++d_id)
			pcl_vols[d_id] = pcl_vols[d_id-1] * 0.25;

		max_depth -= 2;

		double min_pcl_area = pcl_vols[0] * min_area_ratio;
		pcl_hx *= 0.5;
		pcl_hy *= 0.5;

		init_elem_pcl_info(mh);

		// split particles to elements
		PclPart pp, pp_tmp;
		Particle pcl_tmp;
		Particle *pcl2;
		Particle *pcl = mh2ps.first();
		//double vol1 = 0.0; // for debug
		while (mh2ps.not_end_yet(pcl))
		{
			//vol1 += pcl->vol;

			pp_tmp.pe = mh.find_in_which_element(pcl->x, pcl->y);
			if (!pp_tmp.pe)
			{
				pcl2 = pcl;
				pcl = mh2ps.next(pcl);
				mh2ps.del_pcl(*pcl2);
				continue;
			}

			pp_tmp.xl = pcl->x - pcl_hx;
			pp_tmp.yl = pcl->y - pcl_hy;
			pp_tmp.xu = pcl->x + pcl_hx;
			pp_tmp.yu = pcl->y + pcl_hy;
			// if this particle completely lies in elements
			if (mh.is_in_triangle(*pp_tmp.pe, pp_tmp.xl, pp_tmp.yl) &&
				mh.is_in_triangle(*pp_tmp.pe, pp_tmp.xu, pp_tmp.yl) &&
				mh.is_in_triangle(*pp_tmp.pe, pp_tmp.xu, pp_tmp.yu) &&
				mh.is_in_triangle(*pp_tmp.pe, pp_tmp.xl, pp_tmp.yu))
			{
				elem_pcl_infos[pp_tmp.pe->id].add_ori_pcl(*pcl);
				pcl = mh2ps.next(pcl);
				continue;
			}

			// begin splitting particle
			cur_elem_id = pp_tmp.pe->id;
			pp_tmp.xc = pcl->x;
			pp_tmp.yc = pcl->y;
			pp_tmp.depth = 0;

			pp_stack.reset_optimize();
			pp_stack.push(pp_tmp);
			
			pcl2 = pcl;
			pcl = mh2ps.next(pcl);
			mh2ps.del_pcl(*pcl2);

			while (pp_stack.pop(pp))
			{
				if ((mh.is_in_triangle(*pp.pe, pp.xl, pp.yl) &&
					 mh.is_in_triangle(*pp.pe, pp.xu, pp.yl) &&
					 mh.is_in_triangle(*pp.pe, pp.xu, pp.yu) &&
					 mh.is_in_triangle(*pp.pe, pp.xl, pp.yu)) ||
					 pp.depth > max_depth)
				{
					pcl_tmp.x = pp.xc;
					pcl_tmp.y = pp.yc;
					pcl_tmp.vol = pcl_vols[pp.depth];
					add_pcl_prt_to_elem(pp.pe->id, pcl_tmp);
					continue;
				}

				pp_tmp.depth = pp.depth + 1;
				// corner 1
				pp_tmp.xl = pp.xl;
				pp_tmp.xu = pp.xc;
				pp_tmp.yl = pp.yl;
				pp_tmp.yu = pp.yc;
				pp_tmp.xc = (pp_tmp.xl + pp_tmp.xu) * 0.5;
				pp_tmp.yc = (pp_tmp.yl + pp_tmp.yu) * 0.5;
				pp_tmp.pe = mh.find_in_which_element(pp_tmp.xc, pp_tmp.yc);
				if (pp_tmp.pe) pp_stack.push(pp_tmp);
				// corner 2
				pp_tmp.xl = pp.xc;
				pp_tmp.xu = pp.xu;
				pp_tmp.yl = pp.yl;
				pp_tmp.yu = pp.yc;
				pp_tmp.xc = (pp_tmp.xl + pp_tmp.xu) * 0.5;
				pp_tmp.yc = (pp_tmp.yl + pp_tmp.yu) * 0.5;
				pp_tmp.pe = mh.find_in_which_element(pp_tmp.xc, pp_tmp.yc);
				if (pp_tmp.pe) pp_stack.push(pp_tmp);
				// corner 3
				pp_tmp.xl = pp.xc;
				pp_tmp.xu = pp.xu;
				pp_tmp.yl = pp.yc;
				pp_tmp.yu = pp.yu;
				pp_tmp.xc = (pp_tmp.xl + pp_tmp.xu) * 0.5;
				pp_tmp.yc = (pp_tmp.yl + pp_tmp.yu) * 0.5;
				pp_tmp.pe = mh.find_in_which_element(pp_tmp.xc, pp_tmp.yc);
				if (pp_tmp.pe) pp_stack.push(pp_tmp);
				// corner 4
				pp_tmp.xl = pp.xl;
				pp_tmp.xu = pp.xc;
				pp_tmp.yl = pp.yc;
				pp_tmp.yu = pp.yu;
				pp_tmp.xc = (pp_tmp.xl + pp_tmp.xu) * 0.5;
				pp_tmp.yc = (pp_tmp.yl + pp_tmp.yu) * 0.5;
				pp_tmp.pe = mh.find_in_which_element(pp_tmp.xc, pp_tmp.yc);
				if (pp_tmp.pe) pp_stack.push(pp_tmp);
			} // end of splitting particles
			
			// collect particles splitting infos
			merge_pcl_prt(min_pcl_area);
		}

		//double vol2 = 0.0;
		//for (size_t e_id = 0; e_id < elem_num; ++e_id)
		//{
		//	ElemPclInfo &e = elem_pcl_infos[e_id];
		//	for (pcl = e.spt_pcls; pcl; pcl = pcl->next2)
		//		vol2 += pcl->vol;
		//	for (pcl = e.ori_pcls; pcl; pcl = pcl->next2)
		//		vol2 += pcl->vol;
		//}

		merge_particle(min_len_ratio);
		
		//double vol3 = 0.0;
		//for (pcl = mh2ps.first(); mh2ps.not_end_yet(pcl); pcl = mh2ps.next(pcl))
		//	vol3 += pcl->vol;

		//std::cout << "vol1 " << vol1 << "\n";
		//std::cout << "vol2 " << vol2 << "\n";
		//std::cout << "vol3 " << vol3 << "\n";

		clear_elem_pcl_info();
		delete[] pcl_vols;
	}
};

#endif