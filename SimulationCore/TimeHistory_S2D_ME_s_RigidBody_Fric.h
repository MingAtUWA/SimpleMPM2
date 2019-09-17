#ifndef __TIME_HISTORY_S2D_ME_S_RIGIDBODY_FRIC_H__
#define __TIME_HISTORY_S2D_ME_S_RIGIDBODY_FRIC_H__

#include "TimeHistory.h"
#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

/* ===========================================================
Class TimeHistory_S2D_ME_s_RigidBody_Fric
=========================================================== */
class TimeHistory_S2D_ME_s_RigidBody_Fric : public TimeHistory
{
public:
	TimeHistory_S2D_ME_s_RigidBody_Fric() : TimeHistory("S2D_ME_RigidBody_Fric") {}
	~TimeHistory_S2D_ME_s_RigidBody_Fric() {}

	// Output funtion
	int output(void) override { return (*output_func)(*this, *res_file); }

	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = *const_cast<OutputFunc *>(&_res_file.out_func_th_MPM_RigidBody_Fric_PlainBin);
	}

	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = *const_cast<OutputFunc *>(&_res_file.out_func_th_MPM_RigidBody_Fric_XML);
	}
};

#endif