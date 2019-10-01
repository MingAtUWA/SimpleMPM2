#ifndef __RESULT_FILE_PLAIN_BIN_H__
#define __RESULT_FILE_PLAIN_BIN_H__

#include "SimulationCore_pcp.h"

#include <fstream>
#include "ResultFile.h"

// Model: Model_S2D_ME_s_RigidBody
// Step: Step_S2D_ME_s_RigidBody
// TimeHistory: TimeHistory_S2D_ME_s_RigidBody
// ResultFile: ResultFile_PlainBin
class Model_S2D_ME_s_RigidBody;
int rf_out_func_imp_th_MPM_RigidBody_PlainBin(TimeHistory &_th, ResultFile &_rf);

// Model: Model_S2D_ME_s_RigidBody_Fric
// Step: Step_S2D_ME_s_RigidBody_Fric
// TimeHistory: TimeHistory_S2D_ME_s_RigidBody_Fric
// ResultFile: ResultFile_PlainBin
class Model_S2D_ME_s_RigidBody_Fric;
int rf_out_func_imp_th_MPM_RigidBody_Fric_PlainBin(TimeHistory &_th, ResultFile &_rf);

// Model: Model_S2D_CHM_s
// Step: Step_S2D_CHM_s
// TimeHistory: TimeHistory_S2D_CHM_s
// ResultFile: ResultFile_PlainBin
class Model_S2D_CHM_s;
int rf_out_func_imp_th_MPM_CHM_s_PlainBin(TimeHistory &_th, ResultFile &_rf);

class ResultFile_PlainBin : public ResultFile
{
protected:
	std::fstream file;

public:
	ResultFile_PlainBin() {}
	~ResultFile_PlainBin() { finalize(); }
	int init(const char *file_name)
	{
		file.open(file_name, std::ios::binary | std::ios::out);
		if (file.is_open())
			return 0;
		return -1;
	}
	void finalize(void) { file.close(); }

	std::fstream &get_file(void) noexcept { return file; }
//
//public: // functions output model data and time history
//	int output(Model_S2D_ME_s_RigidBody &model);
//	friend int rf_out_func_imp_th_MPM_RigidBody_PlainBin(TimeHistory &_th, ResultFile &_rf);
//	static const OutputFunc out_func_th_MPM_RigidBody_PlainBin;
//
//	int output(Model_S2D_ME_s_RigidBody_Fric &model);
//	friend int rf_out_func_imp_th_MPM_RigidBody_Fric_PlainBin(TimeHistory &_th, ResultFile &_rf);
//	static const OutputFunc out_func_th_MPM_RigidBody_Fric_PlainBin;
//
//	int output(Model_S2D_CHM_s &model);
//	friend int rf_out_func_imp_th_MPM_CHM_s_PlainBin(TimeHistory &_th, ResultFile &_rf);
//	static const OutputFunc out_func_th_MPM_CHM_s_PlainBin;
};

#endif