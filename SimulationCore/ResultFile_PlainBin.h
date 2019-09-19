#ifndef __RESULT_FILE_PLAIN_BIN_H__
#define __RESULT_FILE_PLAIN_BIN_H__

#include "SimulationCore_pcp.h"

#include <fstream>
#include "ResultFile.h"

class Model_S2D_ME_s_RigidBody;
class Model_S2D_ME_s_RigidBody_Fric;

// Model: Model_S2D_ME_s_RigidBody
// Step: Step_S2D_ME_s_RigidBody
// TimeHistory: TimeHistory_S2D_ME_s_RigidBody
// ResultFile: ResultFile_PlainBin
int rf_out_func_imp_th_MPM_RigidBody_PlainBin(TimeHistory &_th, ResultFile &_rf);

// Model: Model_S2D_ME_s_RigidBody_Fric
// Step: Step_S2D_ME_s_RigidBody_Fric
// TimeHistory: TimeHistory_S2D_ME_s_RigidBody_Fric
// ResultFile: ResultFile_PlainBin
int rf_out_func_imp_th_MPM_RigidBody_Fric_PlainBin(TimeHistory &_th, ResultFile &_rf);

class ResultFile_PlainBin : public ResultFile
{
protected:
	std::fstream file;

public:
	ResultFile_PlainBin();
	~ResultFile_PlainBin();
	int init(const char *file_name);
	void finalize(void);

public: // functions output time history
	friend int rf_out_func_imp_th_MPM_RigidBody_PlainBin(TimeHistory &_th, ResultFile &_rf);
	static const OutputFunc out_func_th_MPM_RigidBody_PlainBin;

	friend int rf_out_func_imp_th_MPM_RigidBody_Fric_PlainBin(TimeHistory &_th, ResultFile &_rf);
	static const OutputFunc out_func_th_MPM_RigidBody_Fric_PlainBin;

public: // functions output model data
	int output(Model_S2D_ME_s_RigidBody &model);
	int output(Model_S2D_ME_s_RigidBody_Fric &model);
};

#endif