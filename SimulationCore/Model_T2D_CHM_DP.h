#ifndef __Model_T2D_CHM_DP_H__
#define __Model_T2D_CHM_DP_H__

#include "BC.h"
#include "Model.h"
#include "ModelContainer.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"
#include "DispConRigidCircle.h"

#define N_tol (1.0e-30)

class ConstitutiveModel;

int solve_substep_T2D_CHM_DP(void *_self);

// MPM uses double set of particles
struct Model_T2D_CHM_DP : public Model
{
	friend int solve_substep_T2D_CHM_DP(void *_self);

public: // Node, Element and Particle data structures
	struct Node
	{
		size_t id;
		double x, y;

		// solid phase
		bool has_mp_s;
		double m_s;
		double ax_s, ay_s;
		double vx_s, vy_s;
		double dux_s, duy_s;
		double fx_ext_s, fy_ext_s;
		double fx_int_s, fy_int_s;

		// fluid phase
		bool has_mp_f;
		double m_f;
		double ax_f, ay_f;
		double vx_f, vy_f;
		double dux_f, duy_f;
		double fx_ext_f, fy_ext_f;
		double fx_int_f, fy_int_f;

		// solid - fluid interaction
		double fx_drag, fy_drag;
	
		// has velocity boundary
		// or output
		bool has_vsx_bc;
		bool has_vsy_bc;
		bool has_vfx_bc;
		bool has_vfy_bc;
	};
	
	struct Element;

	struct SolidParticle
	{
		size_t id;
		double x, y;
		double vx, vy;
		double n; // porosity
		double m, density;
		double e11, e22, e12; // strain
		double s11, s22, s12; // effective stress

		ConstitutiveModel *cm;
		inline void set_cm(ConstitutiveModel &_cm) { cm = &_cm; cm->ext_data = this; }

		 // ===== calculation variables =====
		double x_ori, y_ori, ux, uy, vol;

		// shape function value
		double N1, N2, N3;
		
		Element *pe;
		SolidParticle *next_by_elem; // Used by Element
	};

	struct FluidParticle
	{
		size_t id;
		double x, y;
		double vx, vy;
		double m, density;
		double p; // pressure

		ConstitutiveModel *cm;
		inline void set_cm(ConstitutiveModel &_cm) { cm = &_cm; cm->ext_data = this; }

		// ===== calculation variables =====
		bool in_solid; // whehter this fluid particle is in solid
		double x_ori, y_ori, ux, uy, vol;

		// shape function value
		double N1, N2, N3;

		Element *pe;
		FluidParticle *next_by_elem; // Used by Element
	};

	struct Element
	{
		// index
		size_t id;
		//topology
		size_t n1, n2, n3;
		double area, area_2; // 2 * A

		// shape function derivatives at element centre
		double dN1_dx, dN1_dy;
		double dN2_dx, dN2_dy;
		double dN3_dx, dN3_dy;

		// particles list
		SolidParticle *spcls;
		inline void add_pcl(SolidParticle &pcl)
		{
			pcl.next_by_elem = spcls;
			spcls = &pcl;
		}

		FluidParticle *fpcls;
		inline void add_pcl(FluidParticle &pcl)
		{
			pcl.next_by_elem = fpcls;
			fpcls = &pcl;
		}

		// mixed integration
		bool has_mp_ps;
		bool has_mp_pf;
		bool has_mp_m;
		double vol_ps;
		double vol_pf;
		double vol_m; // vol_ms == vol_mf
		double n;
		double s11, s22, s12, p;

		// strain enhancement apprach
		//double dde11, dde22, de12;
		//double de_vol_s, de_vol_f;
	};

public:
	// bg mesh
	size_t node_num;
	Node *nodes;
	size_t elem_num;
	Element *elems;

	// solid particles
	size_t spcl_num;
	SolidParticle *spcls;
	// fluid particles
	size_t fpcl_num;
	FluidParticle *fpcls;

	// boundary conditions
	size_t sbfx_num, sbfy_num;
	BodyForce *sbfxs, *sbfys;
	size_t fbfx_num, fbfy_num;
	BodyForce *fbfxs, *fbfys;
	size_t tx_num, ty_num;
	TractionBC_MPM *txs, *tys;
	size_t sax_num, say_num;
	AccelerationBC *saxs, *says;
	size_t fax_num, fay_num;
	AccelerationBC *faxs, *fays;
	size_t svx_num, svy_num;
	VelocityBC *svxs, *svys;
	size_t fvx_num, fvy_num;
	VelocityBC *fvxs, *fvys;

	double Kf;  // Bulk modulus of water
	double k;   // Permeability
	double miu; // Dynamic viscosity
	double lbd; // bulk viscosity
	// constitutive model
	ModelContainer model_container;

public:
	Model_T2D_CHM_DP();
	~Model_T2D_CHM_DP();
	
	void init_mesh(double *node_coords, size_t n_num,
				   size_t *elem_n_ids,  size_t e_num);
	void init_mesh(TriangleMesh &tri_mesh);
	void clear_mesh(void);

	void init_spcls(size_t num, double _n, double _m, double _density, double _E, double _niu, double _k);
	void init_spcls(TriangleMeshToParticles &mh_2_pcl, double _n, double _density, double _E, double _niu, double _k);
	void alloc_spcls(size_t num);
	void clear_spcls(void);

	void init_fpcls(size_t num, double _m, double _density, double _Kf, double _miu, double _lbd = 0.0);
	void init_fpcls(TriangleMeshToParticles &mh_2_pcl, double _density, double _Kf, double _miu, double _lbd = 0.0);
	void alloc_fpcls(size_t num);
	void clear_fpcls(void);

#define INIT_BC_TEMPLATE(name, type)    \
	void init_ ## name ## s(size_t num) \
	{                                   \
		if (name ## s)                  \
		{                               \
			if (name ## _num < num)		\
				delete[] name ## s;		\
			else                        \
			{                           \
				name ## _num = num;	    \
				return;                 \
			}                           \
		}                               \
		name ## s = new type ## [num];  \
		name ## _num = num;             \
	}                                   \
	void clear_ ## name ## s(void)      \
	{                                   \
		if (name ## s)                  \
		{                               \
			delete[] name ## s;         \
			name ## s = nullptr;        \
		}                               \
		name ## _num = 0;               \
	}

	INIT_BC_TEMPLATE(sbfx, BodyForce)
	INIT_BC_TEMPLATE(sbfy, BodyForce)
	INIT_BC_TEMPLATE(fbfx, BodyForce)
	INIT_BC_TEMPLATE(fbfy, BodyForce)
	INIT_BC_TEMPLATE(tx,  TractionBC_MPM)
	INIT_BC_TEMPLATE(ty,  TractionBC_MPM)
	INIT_BC_TEMPLATE(sax, AccelerationBC)
	INIT_BC_TEMPLATE(say, AccelerationBC)
	INIT_BC_TEMPLATE(fax, AccelerationBC)
	INIT_BC_TEMPLATE(fay, AccelerationBC)
	INIT_BC_TEMPLATE(svx, VelocityBC)
	INIT_BC_TEMPLATE(svy, VelocityBC)
	INIT_BC_TEMPLATE(fvx, VelocityBC)
	INIT_BC_TEMPLATE(fvy, VelocityBC)
	
public:
	template <typename Particle>
	inline bool is_in_triangle(Element &e, Particle &p)
	{
		Node &n1 = nodes[e.n1];
		Node &n2 = nodes[e.n2];
		Node &n3 = nodes[e.n3];
		double a = (n2.x - p.x) * (n3.y - p.y) - (n3.x - p.x) * (n2.y - p.y);
		double b = (n3.x - p.x) * (n1.y - p.y) - (n1.x - p.x) * (n3.y - p.y);
		double c = e.area_2 - a - b;
		bool res = 0.0 <= a && a <= e.area_2
				&& 0.0 <= b && b <= e.area_2
				&& 0.0 <= c && c <= e.area_2;
		if (res)
		{
			p.N1 = a / e.area_2;
			p.N2 = b / e.area_2;
			p.N3 = 1.0 - p.N1 - p.N2;
			if (p.N1 < N_tol)
				p.N1 = N_tol;
			if (p.N2 < N_tol)
				p.N2 = N_tol;
			if (p.N3 < N_tol)
				p.N3 = N_tol;
		}
		return res;
	}

	// ==================== background grid ====================
	// to accelerate searching
protected:
	struct PElement
	{
		Element *e;
		PElement *next;
	};

	struct Grid
	{
		size_t x_id, y_id;
		PElement *pelems;
	};
	
	double grid_x_min, grid_x_max, grid_hx;
	double grid_y_min, grid_y_max, grid_hy;
	size_t grid_x_num, grid_y_num, grid_num;
	Grid *bg_grids;
	MemoryUtilities::ItemBuffer<PElement> pe_buffer;
	
	inline void add_elem_to_grid(Grid &g, Element &e)
	{
		PElement *pe = pe_buffer.alloc();
		pe->e = &e;
		pe->next = g.pelems;
		g.pelems = pe;
	}

	void add_elem_to_bg_grid(Element &e);

public:
	template <typename Point>
	Element *find_in_which_element(Point &p);

	int init_bg_mesh(double hx, double hy);
	inline void clear_bg_mesh(void)
	{
		if (bg_grids)
		{
			delete[] bg_grids;
			bg_grids = nullptr;
		}
		grid_num = 0;
		grid_x_num = 0;
		grid_y_num = 0;
	}
};

template <typename Point>
inline Model_T2D_CHM_DP::Element *Model_T2D_CHM_DP::find_in_which_element(Point &p)
{
	if (p.x < grid_x_min || p.x > grid_x_max ||
		p.y < grid_y_min || p.y > grid_y_max)
		return nullptr;
	size_t x_id = size_t((p.x - grid_x_min) / grid_hx);
	size_t y_id = size_t((p.y - grid_y_min) / grid_hy);
	Grid &g = bg_grids[grid_x_num * y_id + x_id];
	PElement *pelem = g.pelems;
	Element *elem;
	while (pelem)
	{
		elem = pelem->e;
		if (is_in_triangle(*elem, p))
			return elem;
		pelem = pelem->next;
	}
	return nullptr;
}

#undef N_tol

#endif