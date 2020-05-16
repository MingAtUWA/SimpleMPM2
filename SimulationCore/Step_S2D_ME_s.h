#ifndef __Step_S2D_ME_s_h__
#define __Step_S2D_ME_s_h__

#include "ItemArray.hpp"

#include "Step.h"
#include "Model_S2D_ME_s.h"

int solve_substep_S2D_ME_s(void *_self);
int solve_substep_S2D_ME_s_VE(void *_self);

class Step_S2D_ME_s : public Step
{
public:
	typedef Model_S2D_ME_s::ShapeFuncValue ShapeFuncValue;
	typedef Model_S2D_ME_s::Node Node;
	typedef Model_S2D_ME_s::Element Element;
	typedef Model_S2D_ME_s::Particle Particle;

protected:
	int init_calculation(void) override;
	friend int solve_substep_S2D_ME_s(void *_self);
	friend int solve_substep_S2D_ME_s_VE(void *_self);
	int finalize_calculation(void) override;
	
public:
	Step_S2D_ME_s();
	~Step_S2D_ME_s();
	inline void use_volume_enhancement(bool _use = true)
	{
		if (_use)
			solve_substep = &solve_substep_S2D_ME_s_VE;
		else
			solve_substep = &solve_substep_S2D_ME_s;
	}
};

#endif