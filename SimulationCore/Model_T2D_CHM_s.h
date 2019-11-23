#ifndef __Model_T2D_CHM_s_H__
#define __Model_T2D_CHM_s_H__

#include "BC.h"
#include "Model.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"

#define SHAPE_FUNC_VALUE_CONTENT   \
	double N1, N2, N3;             \
	double dN1_dx, dN2_dx, dN3_dx; \
	double dN1_dy, dN2_dy, dN3_dy

#define N1(xi, eta) (xi)
#define N2(xi, eta) (eta)
#define N3(xi, eta) (1.0 - (xi) - (eta))
#define dN1_dxi   (1.0)
#define dN2_dxi   (0.0)
#define dN3_dxi  (-1.0)
#define dN1_deta  (0.0)
#define dN2_deta  (1.0)
#define dN3_deta (-1.0)

struct ShapeFuncValue
{
public:
	SHAPE_FUNC_VALUE_CONTENT;

	inline void cal_shape_func_value(double xi, double eta,
									 double x1, double y1,
									 double x2, double y2,
									 double x3, double y3)
	{
		// calculate shape function
		N1 = N1(xi, eta);
		N2 = N2(xi, eta);
		N3 = N3(xi, eta);
		//std::cout << "func: N1 " << N1 << ", N2 " << N2 << ", N3 " << N3 << "\n";
		// calculate derivatives of shape function
		double dx_dxi, dx_deta, dy_dxi, dy_deta;
		double Jdet;
		double dxi_dx, dxi_dy, deta_dx, deta_dy;
		dx_dxi  = x1 * dN1_dxi  + x2 * dN2_dxi  + x3 * dN3_dxi;
		dx_deta = x1 * dN1_deta + x2 * dN2_deta + x3 * dN3_deta;
		dy_dxi  = y1 * dN1_dxi  + y2 * dN2_dxi  + y3 * dN3_dxi;
		dy_deta = y1 * dN1_deta + y2 * dN2_deta + y3 * dN3_deta;
		Jdet = dx_dxi * dy_deta - dx_deta * dy_dxi;
		dxi_dx  =  dy_deta / Jdet;
		dxi_dy  = -dx_deta / Jdet;
		deta_dx = -dy_dxi  / Jdet;
		deta_dy =  dx_dxi  / Jdet;
		dN1_dx = dN1_dxi * dxi_dx + dN1_deta * deta_dx;
		dN1_dy = dN1_dxi * dxi_dy + dN1_deta * deta_dy;
		dN2_dx = dN2_dxi * dxi_dx + dN2_deta * deta_dx;
		dN2_dy = dN2_dxi * dxi_dy + dN2_deta * deta_dy;
		dN3_dx = dN3_dxi * dxi_dx + dN3_deta * deta_dx;
		dN3_dy = dN3_dxi * dxi_dy + dN3_deta * deta_dy;
		//std::cout << "func: dN1_dx " << dN1_dx << ", dN1_dy " << dN1_dy 
		//		<< ", dN2_dx " << dN2_dx << ", dN2_dy " << dN2_dy
		//		<< ", dN3_dx " << dN3_dx << ", dN3_dy " << dN3_dy << "\n";
	}
};

int solve_substep_T2D_CHM_s(void *_self);

struct Model_T2D_CHM_s : public Model
{
	friend int solve_substep_T2D_CHM_s(void *_self);
public: // Node, Element and Particle data structures
	struct Node
	{
		size_t id;
		double x, y;
		// solid phase
		double m_s;
		double ax_s, ay_s;
		double vx_s, vy_s;
		double dux_s, duy_s;
		double fx_ext_s, fy_ext_s;
		double fx_int_s, fy_int_s;
		// fluid phase
		double m_f;
		double ax_f, ay_f;
		double vx_f, vy_f;
		double dux_f, duy_f;
		double fx_ext_f, fy_ext_f;
		double fx_int_f, fy_int_f;
		// solid - fluid interaction
		double fx_drag, fy_drag;
	};
	
	struct Element
	{
		size_t id, n1, n2, n3; // node index, topology
		double area_2; // 2 * A
		// shape function of element centre
		union { struct { SHAPE_FUNC_VALUE_CONTENT; };  ShapeFuncValue sf; };
		// calculation variables
		double vol, n, s11, s22, s12, p;
	};

	struct Particle
	{
		size_t id;
		double x, y;

		double ux_s, uy_s;
		double vx_s, vy_s;

		double ux_f, uy_f;
		double vx_f, vy_f;

		double n; // porosity
		double m_s; // mass of solid phase 
		double density_s, vol_s;
		double density_f;

		// total strain
		double e11, e22, e12;
		// effective stress
		double s11, s22, s12;
		// pore pressure
		double p;

		// Constitutive model
		double E;   // Elastic modulus
		double niu; // Poisson ratio
		double Kf;  // Bulk modulus of water
		double k;   // Permeability
		double miu; // Dynamic viscosity

		// calculation variables
		double x_ori, y_ori;
		double vol, m_f;
		Element *pe;
		// shape function value
		union { struct { SHAPE_FUNC_VALUE_CONTENT; }; ShapeFuncValue sf; };
	};

public:
	size_t node_num;
	Node *nodes;
	size_t elem_num;
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
	Model_T2D_CHM_s();
	~Model_T2D_CHM_s();
	
	void init_mesh(double *node_coords, size_t n_num,
				   size_t *elem_n_ids,  size_t e_num);
	void clear_mesh(void);

	void init_pcls(size_t num, double n, double m_s, double density_s, double density_f,
				  double E, double niu, double Kf, double k, double miu);
	void clear_pcls(void);

	void init_mesh(TriangleMesh &tri_mesh);
	void init_pcls(TriangleMeshToParticles &mh_2_pcl,
		double n, double density_s, double density_f,
		double E, double niu, double Kf, double k, double miu);

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

	INIT_BC_TEMPLATE(bfx, BodyForce)
	INIT_BC_TEMPLATE(bfy, BodyForce)
	INIT_BC_TEMPLATE(tx, TractionBC_MPM)
	INIT_BC_TEMPLATE(ty, TractionBC_MPM)
	INIT_BC_TEMPLATE(asx, AccelerationBC)
	INIT_BC_TEMPLATE(asy, AccelerationBC)
	INIT_BC_TEMPLATE(afx, AccelerationBC)
	INIT_BC_TEMPLATE(afy, AccelerationBC)
	INIT_BC_TEMPLATE(vsx, VelocityBC)
	INIT_BC_TEMPLATE(vsy, VelocityBC)
	INIT_BC_TEMPLATE(vfx, VelocityBC)
	INIT_BC_TEMPLATE(vfy, VelocityBC)

public:
	inline bool is_in_triangle(Element &elem, double x, double y)
	{
		Node &n1 = nodes[elem.n1];
		Node &n2 = nodes[elem.n2];
		Node &n3 = nodes[elem.n3];
		double a = (n2.x - n1.x) * (y - n1.y) - (x - n1.x) * (n2.y - n1.y);
		double b = (x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (y - n1.y);
		double c = elem.area_2 - a - b;
		return 0.0 <= a && a <= elem.area_2
			&& 0.0 <= b && b <= elem.area_2
			&& 0.0 <= c && c <= elem.area_2;
	}
	
protected:
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
			//std::cout << "p: (" << p.x << ", " << p.y << ") is in e " << e.id << "\n"
			//		  << a/e.area_2 << ", " << b/e.area_2 << "\n";
			ShapeFuncValue &sf = p.sf;
			sf.N1 = a / e.area_2;
			sf.N2 = b / e.area_2;
			sf.N3 = 1.0 - sf.N1 - sf.N2;
			//sf.dN1_dx = (n2.y - n3.y) / e.area_2;
			//sf.dN1_dy = (n3.x - n2.x) / e.area_2;
			//sf.dN2_dx = (n3.y - n1.y) / e.area_2;
			//sf.dN2_dy = (n1.x - n3.x) / e.area_2;
			//sf.dN3_dx = (n1.y - n2.y) / e.area_2;
			//sf.dN3_dy = (n2.x - n1.x) / e.area_2;
			//std::cout << "outf: N1 " << sf.N1 << ", N2 " << sf.N2 << ", N3 " << sf.N3 << "\n";
			//p.sf.cal_shape_func_value(a / e.area_2, b / e.area_2, n1.x, n1.y, n2.x, n2.y, n3.x, n3.y);
		}
		return res;
	}

	// find particle is in which element (to be accelerated)
	// calculate shape function
	// return true if particle is in mesh and false vise versa.
	inline bool init_pcl_cal_var(Particle &pcl)
	{
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			Element &e = elems[e_id];
			if (is_in_triangle(e, pcl))
			{
				pcl.pe = &e;
				return true;
			}
		}
		pcl.pe = nullptr;
		return false;
	}
};

#undef SHAPE_FUNC_VALUE_CONTENT
#undef N1
#undef N2
#undef N3
#undef dN1_dxi
#undef dN2_dxi
#undef dN3_dxi
#undef dN1_deta
#undef dN2_deta
#undef dN3_deta

#endif