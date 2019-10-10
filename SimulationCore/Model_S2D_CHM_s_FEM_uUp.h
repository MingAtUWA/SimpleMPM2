#ifndef __Model_S2D_CHM_S_FEM_uUp_H__
#define __Model_S2D_CHM_S_FEM_uUp_H__

#include "BC.h"
#include "Model.h"

#define N_LOW(xi)  (1.0 - (xi)) / 2.0
#define N_HIGH(xi) (1.0 + (xi)) / 2.0
#define dN_dxi_LOW(xi) -0.5
#define dN_dxi_HIGH(xi) 0.5

namespace Model_S2D_CHM_s_FEM_uUp_Internal
{
	inline double N1(double xi, double eta) noexcept { return N_LOW(xi)  * N_LOW(eta); }
	inline double N2(double xi, double eta) noexcept { return N_HIGH(xi) * N_LOW(eta); }
	inline double N3(double xi, double eta) noexcept { return N_HIGH(xi) * N_HIGH(eta); }
	inline double N4(double xi, double eta) noexcept { return N_LOW(xi)  * N_HIGH(eta); }
};

struct Model_S2D_CHM_s_FEM_uUp : public Model
{
public: // Node, Element and Particle data structures
	struct Node
	{
		size_t index_x, index_y;
		double x, y, p;
		// solid phase
		double ax_s, ay_s;
		double vx_s, vy_s;
		double ux_s, uy_s;
		// fluid phase
		double ax_f, ay_f;
		double vx_f, vy_f;
		double ux_f, uy_f;
	};

#define SHAPE_FUNC_VALUE_CONTENT           \
	double N1, N2, N3, N4;                 \
	double dN1_dx, dN2_dx, dN3_dx, dN4_dx; \
	double dN1_dy, dN2_dy, dN3_dy, dN4_dy
	struct ShapeFuncValue { SHAPE_FUNC_VALUE_CONTENT; };

	struct Element
	{
		size_t index_x, index_y;
		
		double n; // porosity
		double density_s;
		double density_f;

		// Constitutive model
		double E;   // Elastic modulus
		double niu; // Poisson ratio
		double Kf;  // Bulk modulus of water
		double k;   // Permeability
		double miu; // Dynamic viscosity

		// effective stress
		double s11, s22, s12;
		// pore pressure
		double p;
		// total strain
		double e11, e22, e12;

		// shape function value
		union { struct { SHAPE_FUNC_VALUE_CONTENT; }; ShapeFuncValue sf_var; };
	};
#undef SHAPE_FUNC_VALUE_CONTENT

public:
	double h, x0, xn, y0, yn;
	size_t node_x_num, node_y_num, node_num;
	Node *nodes;
	size_t elem_x_num, elem_y_num, elem_num;
	Element *elems;
	double elem_vol;

	// Force BCs (Naumann BCs)
	size_t bfx_num, bfy_num;
	BodyForce *bfxs, *bfys;
	size_t tx_num, ty_num;
	TractionBC_2DFEM *txs, *tys;
	// Displacement BCs (Dirichlet BCs)
	// solid phase
	size_t usx_num, usy_num;
	DisplacementIncBC *usxs, *usys;
	// fluid phase
	size_t ufx_num, ufy_num;
	DisplacementIncBC *ufxs, *ufys;
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

	inline void get_node_id(Element &e,
		size_t &n1_id, size_t &n2_id, size_t &n3_id, size_t &n4_id) const noexcept
	{
		n1_id = e.index_y * node_x_num + e.index_x;
		n2_id = n1_id + 1;
		n3_id = n2_id + node_x_num;
		n4_id = n3_id - 1;
	}
};

#undef N_LOW
#undef N_HIGH
#undef dN_dxi_LOW
#undef dN_dxi_HIGH

#endif