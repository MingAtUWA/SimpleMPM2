#ifndef __Model_S2D_CHM_S_FEM_uUp_H__
#define __Model_S2D_CHM_S_FEM_uUp_H__

#include "BC.h"
#include "Model.h"

#define N_LOW(xi)  (1.0 - (xi)) / 2.0
#define N_HIGH(xi) (1.0 + (xi)) / 2.0
#define dN_dxi_LOW(xi) -0.5
#define dN_dxi_HIGH(xi) 0.5

struct Model_S2D_CHM_s_FEM_uUp : public Model
{
public:
	struct GaussPoint
	{
		size_t index_x, index_y;
		
		double n; // porosity
		double density_s;
		double density_f;

		double ax_s, ay_s;
		double vx_s, vy_s;
		double ux_s, uy_s;
		double ax_f, ay_f;
		double vx_f, vy_f;
		double ux_f, uy_f;

		// effective stress
		double s11, s22, s12;
		// pore pressure
		double p;
		// total strain
		double e11, e22, e12;
		
		// Constitutive model
		double E;   // Elastic modulus
		double niu; // Poisson ratio
		double Kf;  // Bulk modulus of water
		double k;   // Permeability
		double miu; // Dynamic viscosity
	};

	struct Node
	{
		size_t index_x, index_y;
		// solid phase
		double ux_s, uy_s;
		double vx_s, vy_s;
		double ax_s, ay_s;
		// fluid phase
		double ux_f, uy_f;
		double vx_f, vy_f;
		double ax_f, ay_f;
		double p;
	};
	
	struct Element
	{
		size_t index_x, index_y;
		size_t n1_id, n2_id, n3_id, n4_id;
		union
		{
			struct { GaussPoint gp1, gp2, gp3, gp4; };
			GaussPoint gps[4];
		};
	};

	struct ShapeFuncValue
	{
		double N1, N2, N3, N4;
		double dN1_dx, dN2_dx, dN3_dx, dN4_dx;
		double dN1_dy, dN2_dy, dN3_dy, dN4_dy;
	};

public:
	double h, x0, xn, y0, yn;
	size_t elem_x_num, elem_y_num, elem_num;
	Element *elems;
	size_t node_x_num, node_y_num, node_num;
	Node *nodes;

	size_t dof_num;
	ShapeFuncValue gp1_sf, gp2_sf, gp3_sf, gp4_sf;
	double gp_w; // weight = h*h/4
	double dN_dx_mat1[3][8], dN_dx_mat2[3][8];
	double dN_dx_mat3[3][8], dN_dx_mat4[3][8];

	// Force BCs (Naumann BCs)
	size_t bfx_num, bfy_num;
	BodyForce *bfxs, *bfys;
	size_t tx_num, ty_num;
	TractionBC_2DFEM *txs, *tys;
	// Displacement BCs (Dirichlet BCs)
	// solid phase
	size_t usx_num, usy_num;
	DisplacementBC *usxs, *usys;
	// fluid phase
	size_t ufx_num, ufy_num;
	DisplacementBC *ufxs, *ufys;
	size_t pbc_num;
	PressureBC *pbcs;

public:
	Model_S2D_CHM_s_FEM_uUp();
	~Model_S2D_CHM_s_FEM_uUp();

	void init_mesh(double _h, size_t _elem_x_num, size_t _elem_y_num,
				   double x_start = 0.0, double y_start = 0.0);
	void init_mat_param(double n, double density_s, double density_f,
						double E, double niu, double Kf, double k, double miu);
	void clear_mesh(void);

public:
	inline void cal_shape_func_value(ShapeFuncValue &sf_var, double xi, double eta)
	{
		double Nx_low = N_LOW(xi);
		double Nx_high = N_HIGH(xi);
		double Ny_low = N_LOW(eta);
		double Ny_high = N_HIGH(eta);
		sf_var.N1 = Nx_low  * Ny_low;
		sf_var.N2 = Nx_high * Ny_low;
		sf_var.N3 = Nx_high * Ny_high;
		sf_var.N4 = Nx_low  * Ny_high;
		double dNx_dxi_low = dN_dxi_LOW(xi);
		double dNx_dxi_high = dN_dxi_HIGH(xi);
		double dNy_deta_low = dN_dxi_LOW(eta);
		double dNy_deta_high = dN_dxi_HIGH(eta);
		double dxi_dx = 2.0 / h;
		double &deta_dy = dxi_dx;
		sf_var.dN1_dx = dNx_dxi_low  * Ny_low   * dxi_dx;
		sf_var.dN1_dy = Nx_low  * dNy_deta_low  * deta_dy;
		sf_var.dN2_dx = dNx_dxi_high * Ny_low   * dxi_dx;
		sf_var.dN2_dy = Nx_high * dNy_deta_low  * deta_dy;
		sf_var.dN3_dx = dNx_dxi_high * Ny_high  * dxi_dx;
		sf_var.dN3_dy = Nx_high * dNy_deta_high * deta_dy;
		sf_var.dN4_dx = dNx_dxi_low  * Ny_high  * dxi_dx;
		sf_var.dN4_dy = Nx_low  * dNy_deta_high * deta_dy;
	}

#define p_gp1 -0.5773502692
#define p_gp2  0.5773502692
#define N1(xi, eta) (N_LOW(xi)  * N_LOW(eta))
#define N2(xi, eta) (N_HIGH(xi) * N_LOW(eta))
#define N3(xi, eta) (N_HIGH(xi) * N_HIGH(eta))
#define N4(xi, eta) (N_LOW(xi)  * N_HIGH(eta))
	void cal_traction_bc(TractionBC_2DFEM &tbc, double tbc_nf[4])
	{
		double tmp1, tmp2, len;
		double xi_g1, xi_g2, eta_g1, eta_g2, t_g1, t_g2;
		tmp1 = tbc.xi1 - tbc.xi0;
		tmp2 = tbc.eta1 - tbc.eta0;
		len = h * 0.25 * sqrt(tmp1 * tmp1 + tmp2 * tmp2);
		tmp1 = 0.5 * (tbc.xi0 + tbc.xi1);
		tmp2 = 0.5 * (tbc.xi1 - tbc.xi0);
		xi_g1 = tmp1 + p_gp1 * tmp2;
		xi_g2 = tmp1 + p_gp2 * tmp2;
		tmp1 = 0.5 * (tbc.eta0 + tbc.eta1);
		tmp2 = 0.5 * (tbc.eta1 - tbc.eta0);
		eta_g1 = tmp1 + p_gp1 * tmp2;
		eta_g2 = tmp1 + p_gp2 * tmp2;
		tmp1 = 0.5 * (tbc.t0 + tbc.t1);
		tmp2 = 0.5 * (tbc.t1 - tbc.t0);
		t_g1 = tmp1 + p_gp1 * tmp2;
		t_g2 = tmp1 + p_gp2 * tmp2;
		tbc_nf[0] = (N1(xi_g1, eta_g1) * t_g1 + N1(xi_g2, eta_g2) * t_g2) * len;
		tbc_nf[1] = (N2(xi_g1, eta_g1) * t_g1 + N2(xi_g2, eta_g2) * t_g2) * len;
		tbc_nf[2] = (N3(xi_g1, eta_g1) * t_g1 + N3(xi_g2, eta_g2) * t_g2) * len;
		tbc_nf[3] = (N4(xi_g1, eta_g1) * t_g1 + N4(xi_g2, eta_g2) * t_g2) * len;
	}
#undef p_gp1
#undef p_gp2
#undef N1
#undef N2
#undef N3
#undef N4
	
	enum class DOF : size_t
	{
		usx = 0,
		usy = 1,
		ufx = 2,
		ufy = 3,
		p = 4
	};
	inline size_t n_id_to_dof_id(size_t n_id, DOF dof_type) const { return size_t(dof_type) * node_num + n_id; }
};

#undef N_LOW
#undef N_HIGH
#undef dN_dxi_LOW
#undef dN_dxi_HIGH

#endif