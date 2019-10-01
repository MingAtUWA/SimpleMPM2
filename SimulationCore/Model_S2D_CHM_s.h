#ifndef __Model_S2D_CHM_S_H__
#define __Model_S2D_CHM_S_H__

#include "BC.h"
#include "Model.h"

#include "ItemBuffer.hpp"
#include "ItemArray.hpp"

int solve_substep_S2D_CHM_s_GIMP(void *_self);

struct Model_S2D_CHM_s : public Model
{
	friend int solve_substep_S2D_CHM_s_GIMP(void *_self);
public: // Node, Element and Particle data structures
	struct Node
	{
		size_t index_x, index_y;
		// for soil (mixture) phase
		double m_s;
		double ax_s, ay_s;
		double vx_s, vy_s;
		double dux_s, duy_s;
		double fx_ext_m, fy_ext_m;
		double fx_int_m, fy_int_m;
		double fx_kin_f, fy_kin_f;
		// for fluid phase
		double m_tf;
		double ax_f, ay_f;
		double vx_f, vy_f;
		double dux_f, duy_f;
		double fx_ext_tf, fy_ext_tf;
		double fx_int_tf, fy_int_tf;
		double fx_drag_tf, fy_drag_tf;
	};

#define SHAPE_FUNC_VALUE_CONTENT           \
	double N1, N2, N3, N4;                 \
	double dN1_dx, dN2_dx, dN3_dx, dN4_dx; \
	double dN1_dy, dN2_dy, dN3_dy, dN4_dy
	struct ShapeFuncValue { SHAPE_FUNC_VALUE_CONTENT; };

	struct Element
	{
		size_t index_x, index_y;
		// shape function value
		union { struct { SHAPE_FUNC_VALUE_CONTENT; }; ShapeFuncValue sf_var; };

		double avg_s11, avg_s22, avg_s12;
		double pcl_vol;
	};

	// Particle calculation variables
	struct ParticleCalVar
	{
		double x, y, vol;
		size_t elem_x_id, elem_y_id;
		Element *pe;
		Node *pn1, *pn2, *pn3, *pn4;
		// shape function value
		union { struct { SHAPE_FUNC_VALUE_CONTENT; }; ShapeFuncValue sf_var; };
	};
#undef SHAPE_FUNC_VALUE_CONTENT

	struct Particle
	{
		size_t index;

		double x, y;
		double x_ori, y_ori;
		double ux_s, uy_s;
		double vx_s, vy_s;

		double ux_f, uy_f;
		double vx_f, vy_f;

		double n; // porosity
		double m_s; // mass of solid phase 
		double density_s;
		double density_f;

		// effective stress
		double s11, s22, s12;
		// pore pressure
		double p;
		// total strain
		double e11, e22, e12;

		// The number of elements that this particle covers
		// == 0, the particle is out of mesh
		size_t elem_num;
		ParticleCalVar *vars;
		ParticleCalVar var;
		// Position tag
		// == 00000001b at external edge
		// == 00000010b at internal edge
		// == 00000011b at both external and internal edge
		//unsigned char pos_tag;

		// Constitutive model
		double E;   // Elastic modulus
		double niu; // Poisson ratio
		double Kf;  // Bulk modulus of water
		double k;   // Permeability
		double miu; // Dynamic viscosity
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
	size_t asx_num, asy_num;
	AccelerationBC *asxs, *asys;
	size_t vsx_num, vsy_num;
	VelocityBC *vsxs, *vsys;
	// fluid phase
	size_t afx_num, afy_num;
	AccelerationBC *afxs, *afys;
	size_t vfx_num, vfy_num;
	VelocityBC *vfxs, *vfys;

public:
	Model_S2D_CHM_s();
	~Model_S2D_CHM_s();

	void init_mesh(double _h, size_t _elem_x_num, size_t _elem_y_num,
				   double x_start = 0.0, double y_start = 0.0);
	void clear_mesh(void);

	void init_pcl(size_t num, double n, double m_s, double density_s, double density_f,
				  double E, double niu, double Kf, double k, double miu);
	void clear_pcl(void);

public:
#define N_LOW(xi)  (1.0 - (xi)) / 2.0
#define N_HIGH(xi) (1.0 + (xi)) / 2.0
#define dN_dxi_LOW(xi) -0.5
#define dN_dxi_HIGH(xi) 0.5
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
#undef N_LOW
#undef N_HIGH
#undef dN_dxi_LOW
#undef dN_dxi_HIGH

	inline void cal_pcl_shape_func_value(ParticleCalVar &pcl_var)
	{
		pcl_var.pe  = elems + elem_x_num * pcl_var.elem_y_id + pcl_var.elem_x_id;
		pcl_var.pn1 = nodes + node_x_num * pcl_var.elem_y_id + pcl_var.elem_x_id;
		pcl_var.pn2 = pcl_var.pn1 + 1;
		pcl_var.pn3 = pcl_var.pn2 + node_x_num;
		pcl_var.pn4 = pcl_var.pn3 - 1;
		double xi  = 2.0 * ((pcl_var.x - x0) / h - double(pcl_var.elem_x_id)) - 1.0;
		double eta = 2.0 * ((pcl_var.y - y0) / h - double(pcl_var.elem_y_id)) - 1.0;
		cal_shape_func_value(pcl_var.sf_var, xi, eta);
	}

	// ============================ for standard MPM ===============================
public:
	bool init_pcl_standard(Particle &pcl);

	// ================================ for GIMP ===================================
protected:
	MemoryUtilities::ItemBuffer<ParticleCalVar> pcl_var_mem;
	// for particle covers many elements (> 2 elements in either direction)
	struct PclVarInfo {	double pos, len; size_t elem_id; };
	MemoryUtilities::ItemArray<PclVarInfo, 2, 10> x_var_info_buf, y_var_info_buf;
public:
	bool init_pcl_GIMP(Particle &pcl);
};

#endif