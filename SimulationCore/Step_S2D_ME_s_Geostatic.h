#ifndef __Step_S2D_ME_s_Geostatic_h__
#define __Step_S2D_ME_s_Geostatic_h__

#include "ItemArray.hpp"

#include "Step.h"
#include "Model_S2D_ME_s.h"

int solve_substep_S2D_ME_s_Geostatic(void *_self);
int solve_substep_S2D_ME_s_Geostatic_VE(void *_self);

class Step_S2D_ME_s_Geostatic : public Step
{
public:
	typedef Model_S2D_ME_s::ShapeFuncValue ShapeFuncValue;
	typedef Model_S2D_ME_s::Node Node;
	typedef Model_S2D_ME_s::Element Element;
	typedef Model_S2D_ME_s::Particle Particle;

protected:
	double prev_e_kin;

	int init_calculation(void) override;
	friend int solve_substep_S2D_ME_s_Geostatic(void *_self);
	friend int solve_substep_S2D_ME_s_Geostatic_VE(void *_self);
	int finalize_calculation(void) override;
	
public:
	Step_S2D_ME_s_Geostatic();
	~Step_S2D_ME_s_Geostatic();
	inline void use_volume_enhancement(bool _use = true)
	{
		if (_use)
			solve_substep = &solve_substep_S2D_ME_s_Geostatic_VE;
		else
			solve_substep = &solve_substep_S2D_ME_s_Geostatic;
	}

//public:
//	MemoryUtilities::ItemArray<double, 2, 10> pcl_x_len_buf, pcl_y_len_buf;
//	size_t distribute_pcl_mat_to_elems(Particle &pcl);
};

#endif