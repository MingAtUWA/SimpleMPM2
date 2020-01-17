#ifndef __Modified_Cam_Clay_H__
#define __Modified_Cam_Clay_H__

#include <iostream>
#include <cmath>

// This model uses constant stiffness instead of
// deducing it from recompression line
class ModifiedCamClay
{
protected:
	double niu; // possion ratio
	double kappa; // logrithmic recompression modulus
	double lambda; // logrithmic compression modulus
	double fric_angle, M2;
	double e; // void ratio
	double pc; // pre-consolidation stress
	double N, Gamma, M; // normal and critical state line parameters
	union // stress
	{
		double stress[6];
		struct { double s11, s22, s33, s12, s23, s31; };
	};
	union // stress increment
	{
		double dstress[6];
		struct { double ds11, ds22, ds33, ds12, ds23, ds31; };
	};
	union // elastic strain increment
	{
		double dstrain_e[6];
		struct { double dee11, dee22, dee33, dee12, dee23, dee31; };
	};
	union // plastic strain increment
	{
		double dstrain_p[6];
		struct { double dep11, dep22, dep33, dep12, dep23, dep31; };
	};
	union // Elastic stiffness
	{
		double De_mat[6][6];
		double De_mat_array[36];
	};
	union // Elasto-plastic stiffness
	{
		double Dep_mat[6][6];
		double Dep_mat_array[36];
	};

public:
	ModifiedCamClay() : e(0.0), pc(0.0), niu(0.0), 
		kappa(0.0), lambda(0.0), fric_angle(0.0), M2(0.0)
	{
		for (size_t i = 0; i < 6; ++i)
		{
			stress[i] = 0.0;
			dstress[i] = 0.0;
			dstrain_e[i] = 0.0;
			dstrain_p[i] = 0.0;
		}
	}
	~ModifiedCamClay() {}

	// angle unit is degree
	void set_param_NC(double _niu,
		double _kappa, double _lambda, double _fric_ang,
		double _e, double _s[6])
	{
		niu = _niu;
		kappa = _kappa;
		lambda = _lambda;
		fric_angle = _fric_ang / 180.0 * 3.14159265359;
		M = 6.0 * sin(fric_angle) / (3.0 - sin(fric_angle));
		M2 = M * M;
		e = _e;
		// stress
		for (size_t i = 0; i < 6; ++i)
			stress[i] = _s[i];
		// pre-consolidation stress
		double p = cal_p();
		double q = cal_q();
		pc = cal_pc(p, q);
		// stiffness mat
		form_De_mat(p);
		form_Dep_mat();
		// NCL and CSL
		N = e + lambda * log(pc);
		Gamma = N - (lambda - kappa) * log(2.0);
	}

	void set_param_OC(double _niu,
		double _kappa, double _lambda, double _fric_ang,
		double _e, double _s[6], double _pc)
	{
		niu = _niu;
		kappa = _kappa;
		lambda = _lambda;
		fric_angle = _fric_ang / 180.0 * 3.14159265359;
		M2 = 6.0 * sin(fric_angle) / (3.0 - sin(fric_angle));
		M2 = M2 * M2;
		e = _e;
		// stress
		for (size_t i = 0; i < 6; ++i)
			stress[i] = _s[i];
		// pre-consolidation stress
		double p = cal_p();
		double q = cal_q();
		pc = cal_pc(p, q);
		if (_pc > pc) // overconsolidated
		{
			pc = _pc;
			// NCL
			N = e + kappa * log(p) + (lambda - kappa) * log(pc);
			// stiffness mat
			double dg_ds[6], dg_dpc;
			double A_pc, divider;
			form_De_mat(p);
			cal_dg_stress(p, q, dg_ds, dg_dpc);
			A_pc = cal_A_pc(dg_ds);
			divider = cal_divider(dg_ds, dg_dpc, A_pc);
			form_Dep_mat(dg_ds, divider);
		}
		else // normally consolidated
		{
			// stiffness mat
			form_De_mat(p);
			form_Dep_mat();
			// NCL
			N = e + lambda * log(pc);
		}
		// CSL
		Gamma = N - (lambda - kappa) * log(2.0);
	}

	// return value:
	//     = 0 - elstic
	//     > 0 - plastic
	//     < 0 - convergence failure
	int integrate(double dstrain[6])
	{
		for (size_t i = 0; i < 6; ++i)
		{
			// stress
			dstress[i] = De_mat[i][0] * dstrain[0] + De_mat[i][1] * dstrain[1]
					   + De_mat[i][2] * dstrain[2] + De_mat[i][3] * dstrain[3]
					   + De_mat[i][4] * dstrain[4] + De_mat[i][5] * dstrain[5];
			stress[i] += dstress[i];
			// strain
			dstrain_e[i] = dstrain[i];
			dstrain_p[i] = 0.0;
		}

		double p, q, f;
		p = cal_p();
		q = cal_q();
		f = cal_f(p, q);
		if (f <= 0.0)
		{

			form_De_mat(p);
			form_Dep_mat();
			e += (e + 1.0) * (dstrain[0] + dstrain[1] + dstrain[2]);
			return 0;
		}

		double dg_ds[6], dg_dpc;
		double A_pc, divider;
		double dl, dep_cor[6], ds_cor[6];
		for (int iter_id = 0; iter_id < 10; ++iter_id)
		{
			cal_dg_stress(p, q, dg_ds, dg_dpc);
			A_pc = cal_A_pc(dg_ds);
			divider = cal_divider(dg_ds, dg_dpc, A_pc);

			dl = f / divider;
			// strain correction
			dep_cor[0] = dl * dg_ds[0];
			dep_cor[1] = dl * dg_ds[1];
			dep_cor[2] = dl * dg_ds[2];
			dep_cor[3] = dl * dg_ds[3];
			dep_cor[4] = dl * dg_ds[4];
			dep_cor[5] = dl * dg_ds[5];
			// elastic strain
			dstrain_e[0] -= dep_cor[0];
			dstrain_e[1] -= dep_cor[1];
			dstrain_e[2] -= dep_cor[2];
			dstrain_e[3] -= dep_cor[3];
			dstrain_e[4] -= dep_cor[4];
			dstrain_e[5] -= dep_cor[5];
			// plastic strain
			dstrain_p[0] += dep_cor[0];
			dstrain_p[1] += dep_cor[1];
			dstrain_p[2] += dep_cor[2];
			dstrain_p[3] += dep_cor[3];
			dstrain_p[4] += dep_cor[4];
			dstrain_p[5] += dep_cor[5];
			// correct stress
			for (size_t i = 0; i < 6; ++i)
				ds_cor[i] = De_mat[i][0] * dep_cor[0] + De_mat[i][1] * dep_cor[1]
						  + De_mat[i][2] * dep_cor[2] + De_mat[i][3] * dep_cor[3]
						  + De_mat[i][4] * dep_cor[4] + De_mat[i][5] * dep_cor[5];
			dstress[0] -= ds_cor[0];
			dstress[1] -= ds_cor[1];
			dstress[2] -= ds_cor[2];
			dstress[3] -= ds_cor[3];
			dstress[4] -= ds_cor[4];
			dstress[5] -= ds_cor[5];
			stress[0]  -= ds_cor[0];
			stress[1]  -= ds_cor[1];
			stress[2]  -= ds_cor[2];
			stress[3]  -= ds_cor[3];
			stress[4]  -= ds_cor[4];
			stress[5]  -= ds_cor[5];
			// update hardening parameter
			pc += A_pc * dl;

			p = cal_p();
			q = cal_q();
			f = cal_f(p, q);
			if (cal_norm_f(f, p) < 1.0e-8)
			//if (f < 1.0e-4)
			{
				// update e
				e += (e + 1.0) * (dstrain[0] + dstrain[1] + dstrain[2]);
				// new stiffness matrix
				form_De_mat(p);
				cal_dg_stress(p, q, dg_ds, dg_dpc);
				A_pc = cal_A_pc(dg_ds);
				divider = cal_divider(dg_ds, dg_dpc, A_pc);
				form_Dep_mat(dg_ds, divider);
				return iter_id + 1;
			}
		}

		// update e
		e += (e + 1.0) * (dstrain[0] + dstrain[1] + dstrain[2]);
		// new stiffness matrix
		form_De_mat(p);
		cal_dg_stress(p, q, dg_ds, dg_dpc);
		A_pc = cal_A_pc(dg_ds);
		divider = cal_divider(dg_ds, dg_dpc, A_pc);
		form_Dep_mat(dg_ds, divider);
		return -1;
	}

	inline const double *get_stress(void)    noexcept { return stress; }
	inline const double *get_dstress(void)   noexcept { return dstress; }
	inline const double *get_dstrain_e(void) noexcept { return dstrain_e; }
	inline const double *get_dstrain_p(void) noexcept { return dstrain_p; }
	inline const double *get_De_mat(void)  noexcept { return De_mat_array; }
	inline const double *get_Dep_mat(void) noexcept { return Dep_mat_array; }
	inline double get_p(void)  noexcept { return cal_p(); }
	inline double get_q(void)  noexcept { return cal_q(); }
	inline double get_pc(void) noexcept { return pc; }
	inline double get_e_by_strain(void)  noexcept { return e; }
	inline double get_e_by_model(void) noexcept
	{
		return N - lambda * log(pc) + kappa * log(-pc/cal_p());
	}
	inline double get_f(void) noexcept { return cal_f(cal_p(), cal_q()); }
	inline double get_norm_f(void) noexcept { return cal_norm_f(cal_p(), cal_q()); }

protected:
	// form elastic stiffness matrix
	inline void form_De_mat(double p)
	{
		double K_mod = (1.0 + e) * -p / kappa;
		double E_mod = 3.0 * (1.0 - 2.0 * niu) * K_mod;

		double coef;
		coef = (1.0 - niu) / ((1.0 + niu) * (1.0 - 2.0 * niu)) * E_mod;
		De_mat[0][0] = coef;
		De_mat[1][1] = coef;
		De_mat[2][2] = coef;

		coef = niu / ((1.0 + niu) * (1.0 - 2.0 * niu)) * E_mod;
		De_mat[0][1] = coef;
		De_mat[0][2] = coef;
		De_mat[1][0] = coef;
		De_mat[1][2] = coef;
		De_mat[2][0] = coef;
		De_mat[2][1] = coef;
		
		De_mat[0][3] = 0.0;
		De_mat[0][4] = 0.0;
		De_mat[0][5] = 0.0;
		De_mat[1][3] = 0.0;
		De_mat[1][4] = 0.0;
		De_mat[1][5] = 0.0;
		De_mat[2][3] = 0.0;
		De_mat[2][4] = 0.0;
		De_mat[2][5] = 0.0;

		De_mat[3][0] = 0.0;
		De_mat[3][1] = 0.0;
		De_mat[3][2] = 0.0;
		De_mat[4][0] = 0.0;
		De_mat[4][1] = 0.0;
		De_mat[4][2] = 0.0;
		De_mat[5][0] = 0.0;
		De_mat[5][1] = 0.0;
		De_mat[5][2] = 0.0;

		coef = E_mod / (1.0 + niu);
		De_mat[3][3] = coef;
		De_mat[4][4] = coef;
		De_mat[5][5] = coef;

		De_mat[3][4] = 0.0;
		De_mat[3][5] = 0.0;
		De_mat[4][3] = 0.0;
		De_mat[4][5] = 0.0;
		De_mat[5][3] = 0.0;
		De_mat[5][4] = 0.0;
	}
	inline void form_Dep_mat(void) noexcept
	{
		for (size_t i = 0; i < 6; ++i)
			for (size_t j = 0; j < 6; ++j)
				Dep_mat[i][j] = De_mat[i][j];
	}
	inline void form_Dep_mat(double dg_ds[6], double divider) noexcept
	{
		double De_dg_ds[6];
		for (size_t i = 0; i < 6; ++i)
			De_dg_ds[i] = De_mat[i][0] * dg_ds[0] + De_mat[i][1] * dg_ds[1]
					    + De_mat[i][2] * dg_ds[2] + De_mat[i][3] * dg_ds[3]
					    + De_mat[i][4] * dg_ds[4] + De_mat[i][5] * dg_ds[5];
		for (size_t i = 0; i < 6; ++i)
			for (size_t j = 0; j < 6; ++j)
				Dep_mat[i][j] = De_mat[i][j] - De_dg_ds[i] * De_dg_ds[j] / divider;
	}

	inline double cal_p(void) noexcept
	{
		return (s11 + s22 + s33) / 3.0;
	}
	inline double cal_q(void) noexcept
	{
		double s11_s22_diff = s11 - s22;
		double s22_s33_diff = s22 - s33;
		double s33_s11_diff = s33 - s11;
		double q2 = (s11_s22_diff * s11_s22_diff
				   + s22_s33_diff * s22_s33_diff
				   + s33_s11_diff * s33_s11_diff) * 0.5
				  + (s12 * s12 + s23 * s23 + s31 * s31) * 3.0;
		return sqrt(q2);
	}
	// equivalent pre-consolidation stress
	inline double cal_pc(double p, double q) noexcept
	{
		return -(q * q / (M2 * p) + p);
	}
	inline double cal_f(double p, double q) noexcept
	{
		return q * q + M2 * p * (p + pc);
	}
	// normalized f - covergence craterion
	inline double cal_norm_f(double f, double p) noexcept
	{
		return -f / (M2 * p * pc);
	}
	void cal_dg_stress(double p, double q, double dg_ds[6], double &dg_dpc)
	{
		double dg_dp, dg_dq;
		dg_dp = M2 * (2.0 * p + pc);
		dg_dq = 2.0 * q;
		dg_dpc = M2 * p;
		// dq_ds
		double dq_ds11, dq_ds22, dq_ds33;
		double dq_ds12, dq_ds23, dq_ds31;
		dq_ds11 = 0.5 * (s11 + s11 - s22 - s33) / q;
		dq_ds22 = 0.5 * (s22 + s22 - s33 - s11) / q;
		dq_ds33 = 0.5 * (s33 + s33 - s11 - s22) / q;
		dq_ds12 = 3.0 * s12 / q;
		dq_ds23 = 3.0 * s23 / q;
		dq_ds31 = 3.0 * s31 / q;
		// dg_ds
		dg_ds[0] = dg_dp / 3.0 + dg_dq * dq_ds11; // dg_ds11
		dg_ds[1] = dg_dp / 3.0 + dg_dq * dq_ds22; // dg_ds22
		dg_ds[2] = dg_dp / 3.0 + dg_dq * dq_ds33; // dg_ds33
		dg_ds[3] = dg_dq * dq_ds12;               // dg_ds12
		dg_ds[4] = dg_dq * dq_ds23;               // dg_ds23
		dg_ds[5] = dg_dq * dq_ds31;               // dg_ds31
	}
	inline double cal_A_pc(double dg_ds[6]) noexcept
	{
		return -(1.0 + e) * pc / (lambda - kappa) * (dg_ds[0] + dg_ds[1] + dg_ds[2]);
	}
	inline double cal_divider(double dg_ds[6], double dg_dpc, double A_pc) noexcept
	{
		double divider = 0.0;
		for (size_t i = 0; i < 6; i++)
			for (size_t j = 0; j < 6; j++)
				divider += dg_ds[i] * De_mat[i][j] * dg_ds[j];
		divider -= dg_dpc * A_pc;
		return divider;
	}
};

#endif