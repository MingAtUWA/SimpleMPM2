#include "TestsWithGL_pcp.h"

#include <fstream>
#include "ModifiedCamClay.h"
#include "test_sim_core.h"

void test_triaxial_test_cam_clay()
{
	double de = -0.5;
	size_t inc_num = 10000;

	std::fstream res_file;

	const double(*Dep_mat_mcc)[6];
	double Dep_mat[6][6];
	double dstrain[6];
	double Kw_div_n;

	double ini_stress[6] = { -20000.0, -12000.0, -12000.0, 0.0, 0.0, 0.0 };
	double pressure = 0.0;
	
	double Kw = 1.0e6;
	double Kws = 10.0;
	double n = 0.25;
	double pressure_lim = -1000.0;
	ModifiedCamClay mcc;
	//mcc.set_param_NC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress);
	mcc.set_param_OC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress, 50000.0);

	res_file.open("mcc_undrained_res.csv", std::ios::out | std::ios::binary);
	res_file << "strain, p, q, pc, e, e_NCRC, f, res, pressure\n"
		<< 0.0 << ", "
		<< mcc.get_p() << ", "
		<< mcc.get_q() << ", "
		<< mcc.get_pc() << ", "
		<< mcc.get_e_by_strain() << ", "
		<< mcc.get_e_by_model() << ", "
		<< mcc.get_f() << ", 0, "
		<< pressure << "\n";

	de /= double(inc_num);
	for (size_t i = 0; i < inc_num; ++i)
	{
		Dep_mat_mcc = reinterpret_cast<const double(*)[6]>(mcc.get_Dep_mat());
		
		memcpy(Dep_mat, Dep_mat_mcc, sizeof(double)*6*6);
		if (pressure > pressure_lim)
			Kw_div_n = Kw / n;
		else
			Kw_div_n = Kws / n;
		Dep_mat[0][0] -= Kw_div_n;
		Dep_mat[0][1] -= Kw_div_n;
		Dep_mat[0][2] -= Kw_div_n;
		Dep_mat[1][0] -= Kw_div_n;
		Dep_mat[1][1] -= Kw_div_n;
		Dep_mat[1][2] -= Kw_div_n;
		Dep_mat[2][0] -= Kw_div_n;
		Dep_mat[2][1] -= Kw_div_n;
		Dep_mat[2][2] -= Kw_div_n;

		dstrain[0] = de;
		dstrain[1] = -Dep_mat[1][0] / (Dep_mat[1][1] + Dep_mat[1][2]) * de;
		dstrain[2] = -Dep_mat[2][0] / (Dep_mat[2][1] + Dep_mat[2][2]) * de;
		dstrain[3] = 0.0;
		dstrain[4] = 0.0;
		dstrain[5] = 0.0;

		int res = mcc.integrate(dstrain);
		pressure += Kw_div_n * (dstrain[0] + dstrain[1] + dstrain[2]);

		// output result
		res_file << de * (i + 1) << ", "
			<< mcc.get_p() << ", "
			<< mcc.get_q() << ", "
			<< mcc.get_pc() << ", "
			<< mcc.get_e_by_strain() << ", "
			<< mcc.get_e_by_model() << ", "
			<< mcc.get_f() << ", "
			<< res << ", "
			<< pressure << "\n";
	}

	res_file.close();
}
