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

class Model_S2D_ME_s_RigidBody_Fric;
// Model: Model_S2D_ME_s_RigidBody_Fric
// Step: Step_S2D_ME_s_RigidBody_Fric
// TimeHistory: TimeHistory_S2D_ME_s_RigidBody_Fric
// ResultFile: ResultFile_XML
int rf_out_func_imp_th_MPM_RigidBody_Fric_XML(TimeHistory &_th, ResultFile &_rf);

class Model_S2D_CHM_s;
// Model: Model_S2D_CHM_s
// Step: Step_S2D_CHM_s
// TimeHistory: TimeHistory_S2D_CHM_s
// Result: ResultFile_XML
int rf_out_func_imp_th_MPM_CHM_s(TimeHistory &_th, ResultFile &_rf);

class ResultFile_XML : public ResultFile
{
protected:
	std::fstream file;

public:
	ResultFile_XML() {}
	~ResultFile_XML() { finalize(); }
	int init(const char *file_name)
	{
		file.open(file_name, std::ios::binary | std::ios::out);
		if (!file.is_open())
			return -1;

		// xml version info
		const char *file_header = "<?xml version=\"1.0\" encoding=\"ascii\"?>\n"
			"<ResultFile>\n";
		file.write(file_header, strlen(file_header));

		return 0;
	}
	void finalize(void)
	{
		const char *file_ending = "</ResultFile>\n";
		file.write(file_ending, strlen(file_ending));
		file.close();
	}

	std::fstream &get_file(void) noexcept { return file; }

//public: // functions output model
//	int output(Model_S2D_ME_s_RigidBody &model);
//	int output(Model_S2D_ME_s_RigidBody_Fric &model);
//	int output(Model_S2D_CHM_s &model);
//
//public: // functions output time history
//	friend int rf_out_func_imp_th_MPM_RigidBody_XML(TimeHistory &_th, ResultFile &_rf);
//	static const OutputFunc out_func_th_MPM_RigidBody_XML;
//
//	friend int rf_out_func_imp_th_MPM_RigidBody_Fric_XML(TimeHistory &_th, ResultFile &_rf);
//	static const OutputFunc out_func_th_MPM_RigidBody_Fric_XML;
//
//	friend int rf_out_func_imp_th_MPM_CHM_s(TimeHistory &_th, ResultFile &_rf);
//	static const OutputFunc out_func_th_MPM_CHM_s_XML;
};

#endif