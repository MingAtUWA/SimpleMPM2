#ifndef __RESULT_FILE_XML_H__
#define __RESULT_FILE_XML_H__

#include "SimulationCore_pcp.h"

#include <fstream>
#include "ResultFile.h"

class Model_S2D_ME_s_RigidBody;

// Model: Model_S2D_ME_s_RigidBody
// Step: Step_S2D_ME_s_RigidBody
// TimeHistory: TimeHistory_S2D_ME_s_RigidBody
// ResultFile: ResultFile_XML
int rf_out_func_imp_th_MPM_RigidBody_XML(TimeHistory &_th, ResultFile &_rf);

// Model: Model_S2D_ME_s_RigidBody_Fric
// Step: Step_S2D_ME_s_RigidBody_Fric
// TimeHistory: TimeHistory_S2D_ME_s_RigidBody_Fric
// ResultFile: ResultFile_XML
int rf_out_func_imp_th_MPM_RigidBody_Fric_XML(TimeHistory &_th, ResultFile &_rf);

class ResultFile_XML : public ResultFile
{
protected:
	std::fstream file;

public:
	ResultFile_XML();
	~ResultFile_XML();
	int init(const char *file_name);
	void finalize(void);

public: // functions output model
	int output(Model_S2D_ME_s_RigidBody &model);
	int output(Model_S2D_ME_s_RigidBody_Fric &model);

public: // functions output time history
	friend int rf_out_func_imp_th_MPM_RigidBody_XML(TimeHistory &_th, ResultFile &_rf);
	static const OutputFunc out_func_th_MPM_RigidBody_XML;

	friend int rf_out_func_imp_th_MPM_RigidBody_Fric_XML(TimeHistory &_th, ResultFile &_rf);
	static const OutputFunc out_func_th_MPM_RigidBody_Fric_XML;
};

#endif