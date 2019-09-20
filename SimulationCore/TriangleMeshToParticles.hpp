#ifndef __TRIANGLE_MESH_TO_PARTICLES_HPP__
#define __TRIANGLE_MESH_TO_PARTICLES_HPP__

#include "ItemArray.hpp"
#include "geometry.h"
#include "TriangleMesh.h"

// Assumption: ParticleType has member x, y.
template <typename ParticleType>
class TriangleMeshToParticles
{
public:
	enum class GeneratorType : unsigned int
	{
		FirstOrderGuassPoint = 0,
		SecondOrderGaussPoint = 1
	};

public:
	TriangleMeshToParticles(TriangleMesh &_mesh, GeneratorType _type = GeneratorType::FirstOrderGuassPoint) : mesh(_mesh)
	{
		if (unsigned int(_type) < generator_num)
			type = _type;
		else
			type = GeneratorType::FirstOrderGuassPoint; // default value
		cur_generator_func = generator_funcs[unsigned int(_type)];
		cur_generator_pcl_num = generator_pcl_num[unsigned int(_type)];
	}
	~TriangleMeshToParticles() {}
	inline int set_generator(GeneratorType _type) noexcept
	{
		if (unsigned int(_type) < generator_num)
		{
			type = _type;
			cur_generator_func = &generator_funcs[unsigned int(_type)];
			cur_generator_pcl_num = generator_pcl_num[unsigned int(_type)];
			return 0;
		}
		return -1;
	}
	inline size_t get_pcl_num_per_elem(void) const noexcept { return generator_pcl_num[unsigned int(type)]; }
	inline size_t get_pcl_num(void) const noexcept { return mesh.get_elem_num() * generator_pcl_num[unsigned int(type)]; }
	// pcls mush be preallocated
	// treat pcl as square and max_pcl_size is the maximum allowable length of pcl
	// max_pcl_size == 0.0 means no restriction on particle size
	void get_pcls(ParticleType *pcls, double max_pcl_size = 0.0)
	{
		size_t elem_num = mesh.get_elem_num();
		TriangleMesh::Element *elems = mesh.get_elems();
		TriangleMesh::Node *nodes = mesh.get_nodes();
		Point elem_n1, elem_n2, elem_n3;
		if (max_pcl_size == 0.0)
		{
			for (size_t elem_id = 0; elem_id < elem_num; ++elem_id)
			{
				TriangleMesh::Element elem = elems[elem_id];
				TriangleMesh::Node &n1 = nodes[elem.n1];
				elem_n1.x = n1.x;
				elem_n1.y = n1.y;
				TriangleMesh::Node &n2 = nodes[elem.n2];
				elem_n2.x = n2.x;
				elem_n2.y = n2.y;
				TriangleMesh::Node &n3 = nodes[elem.n3];
				elem_n3.x = n3.x;
				elem_n3.y = n3.y;
				(this->*cur_generator_func)(elem_n1, elem_n2, elem_n3, pcls);
				pcls += cur_generator_pcl_num;
			}
		}
		else
		{
			MemoryUtilities::ItemArray<Point, 2, 22> pt_mem;
			Point *low_pts, *upp_pts;
			double dx21, dy21, dx32, dy32;
			for (size_t elem_id = 0; elem_id < elem_num; ++elem_id)
			{
				TriangleMesh::Element elem = elems[elem_id];
				TriangleMesh::Node &n1 = nodes[elem.n1];
				elem_n1.x = n1.x;
				elem_n1.y = n1.y;
				TriangleMesh::Node &n2 = nodes[elem.n2];
				elem_n2.x = n2.x;
				elem_n2.y = n2.y;
				TriangleMesh::Node &n3 = nodes[elem.n3];
				elem_n3.x = n3.x;
				elem_n3.y = n3.y;
				size_t div_num = size_t(floor(sqrt(elem.area)/max_pcl_size));
				if (div_num == 0) div_num = 1;
				dx21 = (elem_n2.x - elem_n1.x) / double(div_num);
				dy21 = (elem_n2.y - elem_n1.y) / double(div_num);
				dx32 = (elem_n3.x - elem_n2.x) / double(div_num);
				dy32 = (elem_n3.y - elem_n2.y) / double(div_num);
				pt_mem.reserve(div_num + div_num + 2);
				upp_pts = pt_mem.get_mem();
				low_pts = upp_pts + div_num + 1;
				low_pts[0].x = elem_n1.x;
				low_pts[0].y = elem_n1.y;
				for (size_t line_id = 0; line_id < div_num; ++line_id)
				{
					Point *tmp = upp_pts;
					upp_pts = low_pts;
					low_pts = tmp;
					low_pts[0].x = upp_pts[0].x + dx21;
					low_pts[0].y = upp_pts[0].y + dy21;
					size_t col_num = line_id + 1;
					for (size_t col_id = 0; col_id < col_num; ++col_id)
					{
						low_pts[col_id + 1].x = low_pts[col_id].x + dx32;
						low_pts[col_id + 1].y = low_pts[col_id].y + dy32;
						(this->*cur_generator_func)(upp_pts[col_id], low_pts[col_id], low_pts[col_id+1], pcls);
						pcls += cur_generator_pcl_num;
					}
					for (size_t col_id = 0; col_id < line_id; ++col_id)
					{
						(this->*cur_generator_func)(upp_pts[col_id], low_pts[col_id], upp_pts[col_id+1], pcls);
						pcls += cur_generator_pcl_num;
					}
				}
			}
		}
	}

protected:
	typedef void(TriangleMeshToParticles::*GeneratorFunc)(Point &p1, Point &p2, Point &p3, ParticleType *pcls);
	static const GeneratorFunc generator_funcs[];
	static const size_t generator_pcl_num[];
	static const size_t generator_num;

	void FirstOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, ParticleType *pcls);
	void SecondOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, ParticleType *pcls);

protected: // run time variables
	TriangleMesh &mesh;
	GeneratorType type;
	GeneratorFunc cur_generator_func;
	size_t cur_generator_pcl_num;
};

template <typename ParticleType>
const typename TriangleMeshToParticles<ParticleType>::GeneratorFunc
	TriangleMeshToParticles<ParticleType>::generator_funcs[] = {
	&FirstOrderGaussPointGenerator,
	&SecondOrderGaussPointGenerator
};

template <typename ParticleType>
const size_t TriangleMeshToParticles<ParticleType>::generator_pcl_num[] = {
	1,
	3
};

template <typename ParticleType>
const size_t TriangleMeshToParticles<ParticleType>::generator_num = 
	sizeof(TriangleMeshToParticles<ParticleType>::generator_funcs) /
	sizeof(TriangleMeshToParticles<ParticleType>::generator_funcs[0]);

template <typename ParticleType>
void TriangleMeshToParticles<ParticleType>::FirstOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, ParticleType *pcls)
{
	pcls[0].x = (p1.x + p2.x + p3.x) / 3.0;
	pcls[0].y = (p1.y + p2.y + p3.y) / 3.0;
}

template <typename ParticleType>
void TriangleMeshToParticles<ParticleType>::SecondOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, ParticleType *pcls)
{
	pcls[0].x = 2.0 / 3.0 * p1.x + 1.0 / 6.0 * p2.x + 1.0 / 6.0 * p3.x;
	pcls[0].y = 2.0 / 3.0 * p1.y + 1.0 / 6.0 * p2.y + 1.0 / 6.0 * p3.y;
	pcls[1].x = 1.0 / 6.0 * p1.x + 1.0 / 6.0 * p2.x + 2.0 / 3.0 * p3.x;
	pcls[1].y = 1.0 / 6.0 * p1.y + 1.0 / 6.0 * p2.y + 2.0 / 3.0 * p3.y;
	pcls[2].x = 1.0 / 6.0 * p1.x + 2.0 / 3.0 * p2.x + 1.0 / 6.0 * p3.x;
	pcls[2].y = 1.0 / 6.0 * p1.y + 2.0 / 3.0 * p2.y + 1.0 / 6.0 * p3.y;
}

#endif