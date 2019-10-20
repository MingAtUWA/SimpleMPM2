#ifndef __Model_S2D_ME_s_up_H__
#define __Model_S2D_ME_s_up_H__

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

namespace Model_S2D_ME_s_up_Internal
{
	inline double N1(double xi, double eta) noexcept { return N_LOW(xi)  * N_LOW(eta); }
	inline double N2(double xi, double eta) noexcept { return N_HIGH(xi) * N_LOW(eta); }
	inline double N3(double xi, double eta) noexcept { return N_HIGH(xi) * N_HIGH(eta); }
	inline double N4(double xi, double eta) noexcept { return N_LOW(xi)  * N_HIGH(eta); }
};

struct Model_S2D_ME_s_up : public Model
{
protected:
	struct ShapeFuncValue { SHAPE_FUNC_VALUE_CONTENT; };

public: // Particle, Node and Element
	struct Node;
	struct Element;
	struct Particle
	{
		size_t index;
		double m;
		double density;
		double x, y;
		double ax, ay;
		double vx, vy;
		double ux, uy;

		double s11, s22, s12, p;
		double e11, e22, e12;

		// Constitutive model
		double E;   // Elastic modulus
		double niu; // Poisson ratio
		double K; // Bulk modulus may = E / ()

		double x_ori, y_ori;
		double vol;
		Element *pe;
		size_t n1_id, n2_id, n3_id, n4_id;
		// shape function value
		union { struct { SHAPE_FUNC_VALUE_CONTENT; }; ShapeFuncValue sf_var; };

		// use by element
		Particle *next;
	};

	struct Node
	{
		size_t index_x, index_y;
		size_t g_id;
	};

	struct Element
	{
		size_t index_x, index_y;

		size_t n1_id, n2_id, n3_id, n4_id;
		// particle list
		Particle *pcls;
		inline void add_pcl(Particle &pcl) { pcl.next = pcls; pcls = &pcl; }
	};

public:
	double h, x0, xn, y0, yn, elem_vol;
	size_t node_x_num, node_y_num, node_num;
	Node *nodes;
	size_t elem_x_num, elem_y_num, elem_num;
	Element *elems;

	size_t pcl_num;
	Particle *pcls;

	// Force BCs (Naumann BCs)
	size_t bfx_num, bfy_num;
	BodyForce *bfxs, *bfys;
	size_t tx_num, ty_num;
	TractionBC_MPM *txs, *tys;
	// Displacement BCs (Dirichlet BCs)
	size_t ux_num, uy_num;
	DisplacementBC *uxs, *uys;
	size_t pbc_num;
	PressureBC *pbcs;

public:
	Model_S2D_ME_s_up();
	~Model_S2D_ME_s_up();

	void init_mesh(double _h, size_t _elem_x_num, size_t _elem_y_num,
				   double x_start = 0.0, double y_start = 0.0);
	void clear_mesh(void);

	void init_pcl(size_t num, double m, double density, double E, double niu, double K);
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

	inline bool init_pcl_cal_var(Particle &pcl)
	{
		if (pcl.x < x0 || pcl.x >= xn || pcl.y < y0 || pcl.y >= yn)
		{
			pcl.pe = nullptr;
			return false;
		}
		size_t elem_x_id = size_t((pcl.x - x0) / h);
		size_t elem_y_id = size_t((pcl.y - y0) / h);
		pcl.pe = elems + elem_x_num * elem_y_id + elem_x_id;
		pcl.n1_id = node_x_num * elem_y_id + elem_x_id;
		pcl.n2_id = pcl.n1_id + 1;
		pcl.n3_id = pcl.n2_id + node_x_num;
		pcl.n4_id = pcl.n3_id - 1;
		double xi  = 2.0 * ((pcl.x - x0) / h - double(elem_x_id)) - 1.0;
		double eta = 2.0 * ((pcl.y - y0) / h - double(elem_y_id)) - 1.0;
		cal_shape_func_value(pcl.sf_var, xi, eta);
		return true;
	}
};

#undef SHAPE_FUNC_VALUE_CONTENT
#undef N_LOW
#undef N_HIGH
#undef dN_dxi_LOW
#undef dN_dxi_HIGH

#endif