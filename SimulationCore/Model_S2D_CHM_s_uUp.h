#ifndef __Model_S2D_CHM_s_uUp_H__
#define __Model_S2D_CHM_s_uUp_H__

#include "BC.h"
#include "Model.h"

#define N_LOW(xi)  (1.0 - (xi)) / 2.0
#define N_HIGH(xi) (1.0 + (xi)) / 2.0
#define dN_dxi_LOW(xi) -0.5
#define dN_dxi_HIGH(xi) 0.5

#define SHAPE_FUNC_VALUE_CONTENT           \
	double N1, N2, N3, N4;                 \
	double dN1_dx, dN2_dx, dN3_dx, dN4_dx; \
	double dN1_dy, dN2_dy, dN3_dy, dN4_dy

namespace Model_S2D_CHM_s_uUp_Internal
{
	inline double N1(double xi, double eta) noexcept { return N_LOW(xi)  * N_LOW(eta); }
	inline double N2(double xi, double eta) noexcept { return N_HIGH(xi) * N_LOW(eta); }
	inline double N3(double xi, double eta) noexcept { return N_HIGH(xi) * N_HIGH(eta); }
	inline double N4(double xi, double eta) noexcept { return N_LOW(xi)  * N_HIGH(eta); }
};

class Step_S2D_CHM_s_uUp;
int solve_substep_S2D_CHM_s_uUp(void *_self);

struct Model_S2D_CHM_s_uUp : public Model
{
	friend Step_S2D_CHM_s_uUp;
	friend int solve_substep_S2D_CHM_s_uUp(void *_self);
public:
	struct ShapeFuncValue { SHAPE_FUNC_VALUE_CONTENT; };
	struct Element;
	struct Particle
	{
		size_t index;

		double x, y;
		double n;
		double m_s;
		double density_s;
		double density_f;
		
		double ax_s, ay_s;
		double vx_s, vy_s;
		double ux_s, uy_s;
		double ax_f, ay_f;
		double ux_f, uy_f;
		double vx_f, vy_f;

		double s11, s22, s12, p;
		double e11, e22, e12;

		double x_ori, y_ori;
		double vol;
		union
		{
			ShapeFuncValue sf;
			struct { SHAPE_FUNC_VALUE_CONTENT; };
		};
		Element *pe;
		Particle *next; // use by Element

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
		size_t g_id;
		double vol, p;
		double ax_s, ay_s, vx_s, vy_s; // soil phase
		double ax_f, ay_f, vx_f, vy_f; // fluid phase
	};

	struct Element
	{
		size_t index_x, index_y;
		size_t n1_id, n2_id, n3_id, n4_id; // node id
		// particle list
		Particle *pcls;
		inline void add_pcl(Particle &pcl) { pcl.next = pcls; pcls = &pcl; }
	};
	
public:
	double h, x0, xn, y0, yn;
	size_t node_x_num, node_y_num, node_num;
	Node *nodes;
	size_t elem_x_num, elem_y_num, elem_num;
	Element *elems;

	size_t pcl_num;
	Particle *pcls;

	// boundary conditions
	size_t bfx_num, bfy_num;
	BodyForce *bfxs, *bfys;
	size_t tx_num, ty_num;
	TractionBC_MPM *txs, *tys;
	// solid phase
	size_t usx_num, usy_num;
	DisplacementBC *usxs, *usys;
	// fluid phase
	size_t ufx_num, ufy_num;
	DisplacementBC *ufxs, *ufys;
	size_t pbc_num;
	PressureBC *pbcs;

public:
	Model_S2D_CHM_s_uUp();
	~Model_S2D_CHM_s_uUp();

	void init_mesh(double _h, size_t _elem_x_num, size_t _elem_y_num,
				   double x_start = 0.0, double y_start = 0.0, double h_tol_ratio = 0.01);
	void clear_mesh(void);

	void init_pcl(size_t num, double n, double m_s, double density_s, double density_f,
				  double E, double niu, double Kf, double k, double miu);
	void clear_pcl(void);

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
		#define deta_dy dxi_dx;
		sf_var.dN1_dx = dNx_dxi_low  * Ny_low   * dxi_dx;
		sf_var.dN1_dy = Nx_low  * dNy_deta_low  * deta_dy;
		sf_var.dN2_dx = dNx_dxi_high * Ny_low   * dxi_dx;
		sf_var.dN2_dy = Nx_high * dNy_deta_low  * deta_dy;
		sf_var.dN3_dx = dNx_dxi_high * Ny_high  * dxi_dx;
		sf_var.dN3_dy = Nx_high * dNy_deta_high * deta_dy;
		sf_var.dN4_dx = dNx_dxi_low  * Ny_high  * dxi_dx;
		sf_var.dN4_dy = Nx_low  * dNy_deta_high * deta_dy;
		#undef deta_dy
	}
	
	inline bool init_pcl_shape_func(Particle &pcl)
	{
		if (pcl.x < x0 || pcl.x >= xn || pcl.y < y0 || pcl.y >= yn)
		{
			pcl.pe = nullptr;
			return false;
		}
		double xi = (pcl.x - x0) / h;
		double eta = (pcl.y - y0) / h;
		size_t elem_x_id = size_t(xi);
		size_t elem_y_id = size_t(eta);
		pcl.pe = elems + elem_x_num * elem_y_id + elem_x_id;
		xi = 2.0 * (xi - double(elem_x_id)) - 1.0;
		eta = 2.0 * (eta - double(elem_y_id)) - 1.0;
		cal_shape_func_value(pcl.sf, xi, eta);
		return true;
	}

public:
	enum class DOF : size_t
	{
		ux_s = 0,
		uy_s = 1,
		ux_f = 2,
		uy_f = 3,
		p = 4
	};
};

#undef SHAPE_FUNC_VALUE_CONTENT
#undef N_LOW
#undef N_HIGH
#undef dN_dxi_LOW
#undef dN_dxi_HIGH

#endif