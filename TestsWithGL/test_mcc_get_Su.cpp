#include "TestsWithGL_pcp.h"

#include <fstream>

#include "ModifiedCamClay.h"

void test_mcc_get_Su(void)
{
	enum class AnalysisType : unsigned char
	{
		TriaxialDrained = 0,
		TriaxialUndrained = 1,
		Consolidation = 2
	};

	AnalysisType tp = AnalysisType::TriaxialUndrained;
	double de = -0.5;
	size_t inc_num = 50000;

	std::fstream res_file;
	const double(*Dep_mat)[6];
	double dstrain[6];
	//const double *stress;

	//double ini_stress[6] = { -18734.0, -18734.0, -18734.0, 0.0, 0.0, 0.0 };
	//double ini_stress[6] = { -20000.0, -12025.0, -12025.0, 0.0, 0.0, 0.0 };
	double ini_stress[6] = { -40000.0, -24050.0, -24050.0, 0.0, 0.0, 0.0 };
	ModifiedCamClay mcc;
	//mcc.set_param_NC(0.3, 0.044, 0.205, 23.5, 2.2516, ini_stress);
	//mcc.set_param_OC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress, 19800.0);
	mcc.set_param_OC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress, 39610.0);
	res_file.open("mcc_Su_res.csv", std::ios::out | std::ios::binary);
	res_file << "strain, p, q, pc, e, e_NCRC, f, res\n"
			 << 0.0 << ", "
			 << mcc.get_p() << ", "
			 << mcc.get_q() << ", "
			 << mcc.get_pc() << ", "
			 << mcc.get_e_by_strain() << ", "
			 << mcc.get_e_by_model() << ", "
			 << mcc.get_f() << ", 0\n";

	de /= double(inc_num);
	for (size_t i = 0; i < inc_num; ++i)
	{
		switch (tp)
		{
		case AnalysisType::TriaxialDrained:
			Dep_mat = reinterpret_cast<const double(*)[6]>(mcc.get_Dep_mat());
			dstrain[0] = de;
			dstrain[1] = -Dep_mat[1][0] / (Dep_mat[1][1] + Dep_mat[1][2]) * de;
			dstrain[2] = -Dep_mat[2][0] / (Dep_mat[2][1] + Dep_mat[2][2]) * de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::TriaxialUndrained:
			dstrain[0] = de;
			dstrain[1] = -0.5 * de;
			dstrain[2] = -0.5 * de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::Consolidation:
			dstrain[0] = de;
			dstrain[1] = 0.0;
			dstrain[2] = 0.0;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		default:
			res_file.close();
			return;
		}

		int res = mcc.integrate(dstrain);

		// output result
		res_file << de * (i + 1) << ", "
				 << mcc.get_p() << ", "
				 << mcc.get_q() << ", "
				 << mcc.get_pc() << ", "
				 << mcc.get_e_by_strain() << ", "
				 << mcc.get_e_by_model() << ", "
				 << mcc.get_f() << ", "
				 << res << "\n";
	}

	res_file.close();
}
