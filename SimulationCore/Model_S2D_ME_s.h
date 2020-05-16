#ifndef __Model_S2D_ME_s_h__
#define __Model_S2D_ME_s_h__

#include "ItemBuffer.hpp"
#include "ItemArray.hpp"

#include "BC.h"
#include "Model.h"

class Step_S2D_ME_s_Geostatic;
int solve_substep_S2D_ME_s_Geostatic(void *_self);
int solve_substep_S2D_ME_s_Geostatic_VE(void *_self); // Volume enhancement

class Step_S2D_ME_s;
int solve_substep_S2D_ME_s(void *_self);
int solve_substep_S2D_ME_s_VE(void *_self);

class Model_S2D_ME_s : public Model
{
	friend class Step_S2D_ME_s_Geostatic;
	friend int solve_substep_S2D_ME_s_Geostatic(void *_self);
	friend int solve_substep_S2D_ME_s_Geostatic_VE(void *_self);

	friend class Step_S2D_ME_s;
	friend int solve_substep_S2D_ME_s(void *_self);
	friend int solve_substep_S2D_ME_s_VE(void *_self);

public:
	struct ShapeFuncValue
	{
		double N1, N2, N3, N4;
		double dN1_dx, dN1_dy;
		double dN2_dx, dN2_dy;
		double dN3_dx, dN3_dy;
		double dN4_dx, dN4_dy;
	};
protected:
	inline void cal_shape_func_value(ShapeFuncValue &sf_var, double xi, double eta)
	{
#define N_LOW(xi)  (1.0 - (xi)) / 2.0
#define N_HIGH(xi) (1.0 + (xi)) / 2.0
#define dN_dxi_LOW -0.5
#define dN_dxi_HIGH 0.5
#define N_tol 1.0e-10
		double Nx_low = N_LOW(xi);
		double Nx_high = N_HIGH(xi);
		double Ny_low = N_LOW(eta);
		double Ny_high = N_HIGH(eta);
		sf_var.N1 = Nx_low  * Ny_low;
		sf_var.N2 = Nx_high * Ny_low;
		sf_var.N3 = Nx_high * Ny_high;
		sf_var.N4 = Nx_low  * Ny_high;
		if (sf_var.N1 < N_tol)
			sf_var.N1 = N_tol;
		if (sf_var.N2 < N_tol)
			sf_var.N2 = N_tol;
		if (sf_var.N3 < N_tol)
			sf_var.N3 = N_tol;
		if (sf_var.N4 < N_tol)
			sf_var.N4 = N_tol;
		double dNx_dxi_low = dN_dxi_LOW;
		double dNx_dxi_high = dN_dxi_HIGH;
		double dNy_deta_low = dN_dxi_LOW;
		double dNy_deta_high = dN_dxi_HIGH;
		double dxi_dx = 2.0 / hx;
		double deta_dy = 2.0 / hy;
		sf_var.dN1_dx = dNx_dxi_low  * Ny_low   * dxi_dx;
		sf_var.dN1_dy = Nx_low  * dNy_deta_low  * deta_dy;
		sf_var.dN2_dx = dNx_dxi_high * Ny_low   * dxi_dx;
		sf_var.dN2_dy = Nx_high * dNy_deta_low  * deta_dy;
		sf_var.dN3_dx = dNx_dxi_high * Ny_high  * dxi_dx;
		sf_var.dN3_dy = Nx_high * dNy_deta_high * deta_dy;
		sf_var.dN4_dx = dNx_dxi_low  * Ny_high  * dxi_dx;
		sf_var.dN4_dy = Nx_low  * dNy_deta_high * deta_dy;
#undef N_LOW
#undef N_HIGH
#undef dN_dxi_LOW
#undef dN_dxi_HIGH
#undef N_tol
	}
	
public:
	struct Node
	{
		size_t x_id, y_id;
		bool has_mp;
		double m;
		double ax, ay;
		double vx, vy;
		double dux, duy;
		double fx_ext, fy_ext;
		double fx_int, fy_int;
		// strain enhancement
		double vol;
		double de_vol;
	};
	
	struct Element
	{
		size_t x_id, y_id;
		Node *pn1, *pn2, *pn3, *pn4;
		// mixed integration
		bool has_pcl;
		double mi_vol;
		double s11, s22, s12;
		// volume enhancement
		double ve_vol;
		// strain enhancement
		double de_vol;
	};

	struct Particle
	{
		size_t id;
		double x, y;
		double vx, vy;
		double m, density;
		double s11, s22, s12;
		double e11, e22, e12;
		double ar; // aspective ratio
		// calculation variables
		double x_ori, y_ori;
		double ux, uy;
		double vol;
		Element *pe;
		ShapeFuncValue sfv;
	};

protected:
	double x0, xn, hx;
	double y0, yn, hy;
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
	size_t ax_num, ay_num;
	AccelerationBC *axs, *ays;
	size_t vx_num, vy_num;
	VelocityBC *vxs, *vys;

	// calculation variables
	double h_tol;
	double elem_area;
	// shape function at element centre
	ShapeFuncValue elem_sfv;
	
	// Constitutive model
	double E;   // Elastic modulus
	double niu; // Poisson ratio
	
public:
	Model_S2D_ME_s();
	~Model_S2D_ME_s();

	int init_mesh(double _x0, double _y0, double _xn, double _yn,
				  size_t _x_num, size_t _y_num, double h_tol_r = 0.001);
	void clear_mesh();

	void alloc_pcls(size_t num);
	void init_pcls(size_t num, double m, double density, double _ar = 1.0);
	void init_pcls(double _x0, double _y0, double _xn, double _yn,
				   size_t _x_num, size_t _y_num, double density);
	void clear_pcls();

	void init_mat_model_le(double _E, double _niu) { E = _E, niu = _niu; }

	inline size_t get_node_num() { return node_num; }
	inline Node *get_nodes() { return nodes; }
	inline size_t get_elem_num() { return elem_num; }
	inline Element *get_elems() { return elems; }
	inline size_t get_pcl_num() { return pcl_num; }
	inline Particle *get_pcls() { return pcls; }

	inline double get_hx() { return hx; }
	inline double get_hy() { return hy; }
	inline double get_x0() { return x0; }
	inline double get_y0() { return y0; }
	inline double get_xn() { return xn; }
	inline double get_yn() { return yn; }
	inline size_t get_elem_x_num() { return elem_x_num; }
	inline size_t get_elem_y_num() { return elem_y_num; }
	inline double get_node_x(Node &n) { return x0 + double(n.x_id) * hx; }
	inline double get_node_y(Node &n) { return y0 + double(n.y_id) * hy; }

#define INIT_BC_TEMPLATE(name, type)    \
	void clear_ ## name ## s()          \
	{                                   \
		if (name ## s)                  \
		{                               \
			delete[] name ## s;         \
			name ## s = nullptr;        \
		}                               \
		name ## _num = 0;               \
	}                                   \
	void init_ ## name ## s(size_t num) \
	{                                   \
		clear_ ## name ## s();          \
		if (name ## s)                  \
		{                               \
			if (name ## _num < num)     \
			{                           \
				delete[] name##s;       \
			}                           \
			else                        \
			{                           \
				name ## _num = num;     \
				return;                 \
			}                           \
		}                               \
		name ## s = new type ## [num];  \
		name ## _num = num;             \
	}                                   \
	inline size_t get_ ## name ## _num() { return name ## _num; } \
	inline type *get_ ## name ## s() { return name ## s; }
	
	INIT_BC_TEMPLATE(bfx, BodyForce)
	INIT_BC_TEMPLATE(bfy, BodyForce)
	INIT_BC_TEMPLATE(tx, TractionBC_MPM)
	INIT_BC_TEMPLATE(ty, TractionBC_MPM)
	INIT_BC_TEMPLATE(ax, AccelerationBC)
	INIT_BC_TEMPLATE(ay, AccelerationBC)
	INIT_BC_TEMPLATE(vx, VelocityBC)
	INIT_BC_TEMPLATE(vy, VelocityBC)
#undef INIT_BC_TEMPLATE

protected: // helper functions
	bool find_in_which_element(Particle &pcl);
public:
	MemoryUtilities::ItemArray<double, 2, 10> pcl_x_len_buf, pcl_y_len_buf;
	size_t distribute_pcl_mat_to_elems(Particle &pcl);
};

inline bool Model_S2D_ME_s::find_in_which_element(Particle &pcl)
{
	// check if particles is out of mesh
	if (pcl.x <= x0 || pcl.x >= xn || pcl.y <= y0 || pcl.y >= yn)
	{
		pcl.pe = nullptr;
		return false;
	}

	size_t elem_x_id, elem_y_id;
	elem_x_id = size_t((pcl.x - x0) / hx);
	elem_y_id = size_t((pcl.y - y0) / hy);
	pcl.pe = elems + elem_x_num * elem_y_id + elem_x_id;

	double xi = 2.0 * ((pcl.x - x0) / hx - double(elem_x_id)) - 1.0;
	double eta = 2.0 * ((pcl.y - y0) / hy - double(elem_y_id)) - 1.0;
	cal_shape_func_value(pcl.sfv, xi, eta);

	return true;
}

#endif