#ifndef __MODEL_S2D_ME_S_RIGIDBODY_FRIC_H__
#define __MODEL_S2D_ME_S_RIGIDBODY_FRIC_H__

#include "ItemArray.hpp"
#include "ContactState.h"
#include "RigidBody.h"
#include "BC.h"

#include "Model.h"

struct Model_S2D_ME_s_RigidBody_Fric : public Model
{
public:
	struct Node;
	struct Element;

	struct Particle
	{
		size_t index;

		double x, y;
		double x_ori, y_ori;
		double ux, uy;

		double vx, vy;
		double m, density;

		double s11, s22, s12;
		double e11, e22, e12;

		ContactState contact_state;

		// constitutive model
		double E, niu;

		// calculation variables
		double vol;

		// location
		Element *pelem;
		Node *pn1, *pn2, *pn3, *pn4;
		double N1, N2, N3, N4;
		double dN1_dx, dN2_dx, dN3_dx, dN4_dx;
		double dN1_dy, dN2_dy, dN3_dy, dN4_dy;
		
		// for contact detection
		Particle *next_by_elem;
	};

	struct Node
	{
		size_t index_x, index_y;
		double m;
		double ax, ay;
		double vx, vy;
		double dux, duy;
		double fx_ext, fy_ext;
		double fx_int, fy_int;
		double fx_con, fy_con; // contact force
	};

	// element is designed for contact detection
	struct Element
	{
		size_t index_x, index_y;
		Particle *pcls;
		bool has_rigid_object;
		inline void reset(void) noexcept
		{
			pcls = nullptr;
			has_rigid_object = false;
		}
		inline void add_pcl(Particle &pcl) noexcept
		{
			pcl.next_by_elem = pcls;
			pcls = &pcl;
		}
	};

public: // Model data
	// background mesh
	double h;
	double x0, xn, y0, yn;
	size_t elem_x_num, elem_y_num, elem_num;
	Element *elems;
	size_t node_x_num, node_y_num, node_num;
	Node *nodes;
	// particles
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

public:
	Model_S2D_ME_s_RigidBody_Fric();
	~Model_S2D_ME_s_RigidBody_Fric();

	void init_mesh(double grid_size, size_t _elem_x_num, size_t _elem_y_num,
				   double x_start = 0.0, double y_start = 0.0);
	void clear_mesh(void);
	void init_pcl(size_t num, double m, double density, double E, double niu);
	void clear_pcl(void);

	size_t get_node_num(void) const noexcept { return node_num; }
	size_t get_pcl_num(void) const noexcept { return pcl_num; }

	inline bool init_pcl_cal_var(Particle &pcl) noexcept
	{
		if (pcl.x < x0 || pcl.x >= xn ||
			pcl.y < y0 || pcl.y >= yn)
		{
			pcl.pelem = nullptr;
			return false;
		}

		// init the element where the particle lies
		size_t e_x_id = size_t((pcl.x - x0) / h);
		size_t e_y_id = size_t((pcl.y - y0) / h);
		pcl.pelem = elems + elem_x_num * e_y_id + e_x_id;
		pcl.pn1 = nodes + node_x_num * e_y_id + e_x_id;
		pcl.pn2 = pcl.pn1 + 1;
		pcl.pn3 = pcl.pn2 + node_x_num;
		pcl.pn4 = pcl.pn3 - 1;
		
		// init shape function
#define N_LOW(xi)  (1.0 - (xi)) / 2.0
#define N_HIGH(xi) (1.0 + (xi)) / 2.0
#define dN_dxi_LOW(xi) -0.5
#define dN_dxi_HIGH(xi) 0.5
		double xi = 2.0 * ((pcl.x - x0) / h - double(e_x_id)) - 1.0;
		double Nx_low = N_LOW(xi);
		double Nx_high = N_HIGH(xi);
		double dNx_dxi_low = dN_dxi_LOW(xi);
		double dNx_dxi_high = dN_dxi_HIGH(xi);
		double eta = 2.0 * ((pcl.y - y0) / h - double(e_y_id)) - 1.0;
		double Ny_low = N_LOW(eta);
		double Ny_high = N_HIGH(eta);
		double dNy_deta_low = dN_dxi_LOW(eta);
		double dNy_deta_high = dN_dxi_HIGH(eta);
		pcl.N1 = Nx_low  * Ny_low;
		pcl.N2 = Nx_high * Ny_low;
		pcl.N3 = Nx_high * Ny_high;
		pcl.N4 = Nx_low  * Ny_high;
		double dxi_dx = 2.0 / h; // = deta_dy
		pcl.dN1_dx = dNx_dxi_low  * Ny_low   * dxi_dx;
		pcl.dN1_dy = Nx_low  * dNy_deta_low  * dxi_dx;
		pcl.dN2_dx = dNx_dxi_high * Ny_low   * dxi_dx;
		pcl.dN2_dy = Nx_high * dNy_deta_low  * dxi_dx;
		pcl.dN3_dx = dNx_dxi_high * Ny_high  * dxi_dx;
		pcl.dN3_dy = Nx_high * dNy_deta_high * dxi_dx;
		pcl.dN4_dx = dNx_dxi_low  * Ny_high  * dxi_dx;
		pcl.dN4_dy = Nx_low  * dNy_deta_high * dxi_dx;
#undef N_LOW
#undef N_HIGH
#undef dN_dxi_LOW
#undef dN_dxi_HIGH

		// init particle volume
		pcl.vol = pcl.m / pcl.density;
		
		return true;
	}

public:
	RigidBody rigid_body;

public: // helper data and functions
	typedef MemoryUtilities::ItemArray<long long, 2, 30> PreAllocIdArray;
	struct LongLongRange { long long lower, upper; };
	LongLongRange y_id_range1, y_id_range2, y_id_range3, y_id_range4;
	PreAllocIdArray x_ids_mem1, x_ids_mem2, x_ids_mem3, x_ids_mem4;
	// get intersection point (y_id, x_id)
	// return false is y_id is not in the range of background mesh
	bool get_intersect_points(double x1, double y1, double x2, double y2,
					LongLongRange &y_id_range, PreAllocIdArray &x_ids_mem);
	// update the range of x_id for each y_id
	long long y_id_min;
	MemoryUtilities::ItemArray<LongLongRange, 2, 30> x_range_mem;
	void update_x_id_range(LongLongRange &y_id_range, PreAllocIdArray & x_ids_mem);
public:
	void rasterize_rect_on_grid(double x1, double y1, double x2, double y2,
								double x3, double y3, double x4, double y4);

public:
	int get_pcls_from_mesh(const char *mesh_file_name,
		double density, double E, double niu, double max_pcl_area = 0.0);
};

#endif