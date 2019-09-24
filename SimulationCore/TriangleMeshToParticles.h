#ifndef __TRIANGLE_MESH_TO_PARTICLES_HPP__
#define __TRIANGLE_MESH_TO_PARTICLES_HPP__

#include "ItemArray.hpp"
#include "ItemBuffer.hpp"

#include "geometry.h"
#include "TriangleMesh.h"

class TriangleMeshToParticles
{
public:
	struct Particle
	{
		double x, y, vol;
	protected:
		friend TriangleMeshToParticles;
		Particle *prev;
	};
	typedef MemoryUtilities::ItemBuffer<Particle> ParticleBuffer;
	enum class GeneratorType : unsigned int
	{
		FirstOrderGaussPoint = 0,
		SecondOrderGaussPoint = 1
	};
protected: // particle generator functions
	typedef void(TriangleMeshToParticles::*GeneratorFunc)(Point &p1, Point &p2, Point &p3, double vol);

public:
	inline Particle *first(void) noexcept { return top; }
	inline Particle *next(Particle *cur) noexcept { return cur->prev; }

public:
	TriangleMeshToParticles(TriangleMesh &_mesh,
		GeneratorType _type = GeneratorType::FirstOrderGaussPoint) : 
		mesh(_mesh), top(nullptr), pcl_num(0)
	{
		if (unsigned int(_type) < generator_num)
			type = _type;
		else
			type = GeneratorType::FirstOrderGaussPoint; // default value
		cur_generator_func = generator_funcs[unsigned int(_type)];
		cur_generator_pcl_num = generator_pcl_num[unsigned int(_type)];
	}
	~TriangleMeshToParticles() { clear(); }
	inline void clear(void) noexcept
	{
		top = nullptr;
		pcl_num = 0;
		particle_buffer.clear();
	}
	inline int set_generator(GeneratorType _type) noexcept
	{
		if (unsigned int(_type) < generator_num)
		{
			type = _type;
			cur_generator_func = generator_funcs[unsigned int(_type)];
			cur_generator_pcl_num = generator_pcl_num[unsigned int(_type)];
			return 0;
		}
		return -1;
	}

public: // main function
	// max_pcl_size == 0.0 means no restriction on particle size
	void generate_pcls(double max_pcl_area = 0.0);
	inline size_t get_pcl_num_per_elem(void) const noexcept { return generator_pcl_num[unsigned int(type)]; }
	inline size_t get_pcl_num(void) const noexcept { return pcl_num; }

protected:
	TriangleMesh &mesh;
	GeneratorType type;
	GeneratorFunc cur_generator_func;
	size_t cur_generator_pcl_num;

	Particle *top;
	size_t pcl_num;
	ParticleBuffer particle_buffer;
	inline Particle *add_pcl(void)
	{
		Particle *res = particle_buffer.alloc();
		res->prev = top;
		top = res;
		++pcl_num;
		return res;
	}

protected:
	static const GeneratorFunc generator_funcs[];
	static const size_t generator_pcl_num[];
	static const size_t generator_num;
	// particle generator functions
	void FirstOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, double vol);
	void SecondOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, double vol);
};

#endif