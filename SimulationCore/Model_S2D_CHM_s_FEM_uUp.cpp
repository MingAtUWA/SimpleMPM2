#include "SimulationCore_pcp.h"

#include "Model_S2D_CHM_s_FEM_uUp.h"

Model_S2D_CHM_s_FEM_uUp::Model_S2D_CHM_s_FEM_uUp() :
	Model("Model_S2D_CHM_s_FEM_uUp"),
	elems(nullptr), elem_x_num(0), elem_y_num(0), elem_num(0),
	nodes(nullptr), node_x_num(0), node_y_num(0), node_num(0),
	bfx_num(0), bfy_num(0), bfxs(nullptr), bfys(nullptr),
	tx_num(0),  ty_num(0),  txs(nullptr),  tys(nullptr),
	usx_num(0), usy_num(0), usxs(nullptr), usys(nullptr),
	ufx_num(0), ufy_num(0), ufxs(nullptr), ufys(nullptr),
	pbc_num(0), pbcs(nullptr) {}

Model_S2D_CHM_s_FEM_uUp::~Model_S2D_CHM_s_FEM_uUp()
{
	clear_mesh();
	if (bfxs)
	{
		delete[] bfxs;
		bfxs = nullptr;
		bfx_num = 0;
	}
	if (bfys)
	{
		delete[] bfys;
		bfys = nullptr;
		bfy_num = 0;
	}
	if (txs)
	{
		delete[] txs;
		txs = nullptr;
		tx_num = 0;
	}
	if (tys)
	{
		delete[] tys;
		tys = nullptr;
		ty_num = 0;
	}
	if (usxs)
	{
		delete[] usxs;
		usxs = nullptr;
		usx_num = 0;
	}
	if (usys)
	{
		delete[] usys;
		usys = nullptr;
		usy_num = 0;
	}
	if (ufxs)
	{
		delete[] ufxs;
		ufxs = nullptr;
		ufx_num = 0;
	}
	if (ufys)
	{
		delete[] ufys;
		ufys = nullptr;
		ufy_num = 0;
	}
	if (pbcs)
	{
		delete[] pbcs;
		pbcs = nullptr;
		pbc_num = 0;
	}
}

void Model_S2D_CHM_s_FEM_uUp::init_mesh(
	double _h, size_t _elem_x_num, size_t _elem_y_num,
	double x_start, double y_start)
{
	clear_mesh();
	h = _h;
	x0 = x_start;
	xn = x0 + _h * double(_elem_x_num);
	y0 = y_start;
	yn = y0 + _h * double(_elem_y_num);

	size_t i, j, k;
	node_x_num = _elem_x_num + 1;
	node_y_num = _elem_y_num + 1;
	node_num = node_x_num * node_y_num;
	nodes = new Node[node_num];
	k = 0;
	for (i = 0; i < node_y_num; ++i)
		for (j = 0; j < node_x_num; ++j)
		{
			Node &n = nodes[k++];
			n.index_x = j;
			n.index_y = i;
		}
	
	elem_x_num = _elem_x_num;
	elem_y_num = _elem_y_num;
	elem_num = elem_x_num * elem_y_num;
	elems = new Element[elem_num];
	k = 0;
	for (i = 0; i < elem_y_num; ++i)
		for (j = 0; j < elem_x_num; ++j)
		{
			Element &e = elems[k++];
			e.index_x = j;
			e.index_y = i;
			e.n1_id = node_x_num * i + j;
			e.n2_id = e.n1_id + 1;
			e.n3_id = e.n2_id + node_x_num;
			e.n4_id = e.n3_id - 1;
		}

	dof_num = node_num * 5;
	cal_shape_func_value(gp1_sf, -0.5773502692, -0.5773502692);
	cal_shape_func_value(gp2_sf,  0.5773502692, -0.5773502692);
	cal_shape_func_value(gp3_sf,  0.5773502692,  0.5773502692);
	cal_shape_func_value(gp4_sf, -0.5773502692,  0.5773502692);
	gp_w = h * h * 0.25;
#define Form_dN_dx_mat(id)                              \
	dN_dx_mat ## id ## [0][0] = gp ## id ## _sf.dN1_dx; \
	dN_dx_mat ## id ## [0][1] = gp ## id ## _sf.dN2_dx; \
	dN_dx_mat ## id ## [0][2] = gp ## id ## _sf.dN3_dx; \
	dN_dx_mat ## id ## [0][3] = gp ## id ## _sf.dN4_dx; \
	dN_dx_mat ## id ## [0][4] = 0.0;                    \
	dN_dx_mat ## id ## [0][5] = 0.0;                    \
	dN_dx_mat ## id ## [0][6] = 0.0;                    \
	dN_dx_mat ## id ## [0][7] = 0.0;                    \
	dN_dx_mat ## id ## [1][0] = 0.0;                    \
	dN_dx_mat ## id ## [1][1] = 0.0;                    \
	dN_dx_mat ## id ## [1][2] = 0.0;                    \
	dN_dx_mat ## id ## [1][3] = 0.0;                    \
	dN_dx_mat ## id ## [1][4] = gp ## id ## _sf.dN1_dy; \
	dN_dx_mat ## id ## [1][5] = gp ## id ## _sf.dN2_dy; \
	dN_dx_mat ## id ## [1][6] = gp ## id ## _sf.dN3_dy; \
	dN_dx_mat ## id ## [1][7] = gp ## id ## _sf.dN4_dy; \
	dN_dx_mat ## id ## [2][0] = gp ## id ## _sf.dN1_dy; \
	dN_dx_mat ## id ## [2][1] = gp ## id ## _sf.dN2_dy; \
	dN_dx_mat ## id ## [2][2] = gp ## id ## _sf.dN3_dy; \
	dN_dx_mat ## id ## [2][3] = gp ## id ## _sf.dN4_dy; \
	dN_dx_mat ## id ## [2][4] = gp ## id ## _sf.dN1_dx; \
	dN_dx_mat ## id ## [2][5] = gp ## id ## _sf.dN2_dx; \
	dN_dx_mat ## id ## [2][6] = gp ## id ## _sf.dN3_dx; \
	dN_dx_mat ## id ## [2][7] = gp ## id ## _sf.dN4_dx
	Form_dN_dx_mat(1);
	Form_dN_dx_mat(2);
	Form_dN_dx_mat(3);
	Form_dN_dx_mat(4);
}

void Model_S2D_CHM_s_FEM_uUp::init_mat_param(
	double n, double density_s, double density_f,
	double E, double niu, double Kf, double k, double miu)
{
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		n.ux_s = 0.0;
		n.uy_s = 0.0;
		n.ux_f = 0.0;
		n.uy_f = 0.0;
		n.p = 0.0;
	}
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &e = elems[e_id];
#define Init_Gauss_Point_Mat_Param(id)        \
		GaussPoint &gp ## id = e.gp ## id ##; \
		gp ## id ##.n = n;      \
		gp ## id ##.density_s = density_s; \
		gp ## id ##.density_f = density_f; \
		gp ## id ##.ax_s = 0.0; \
		gp ## id ##.ay_s = 0.0; \
		gp ## id ##.vx_s = 0.0; \
		gp ## id ##.vy_s = 0.0; \
		gp ## id ##.ux_s = 0.0; \
		gp ## id ##.uy_s = 0.0; \
		gp ## id ##.ax_f = 0.0; \
		gp ## id ##.ay_f = 0.0; \
		gp ## id ##.vx_f = 0.0; \
		gp ## id ##.vy_f = 0.0; \
		gp ## id ##.ux_f = 0.0; \
		gp ## id ##.uy_f = 0.0; \
		gp ## id ##.s11 = 0.0;  \
		gp ## id ##.s22 = 0.0;  \
		gp ## id ##.s12 = 0.0;  \
		gp ## id ##.p = 0.0;    \
		gp ## id ##.e11 = 0.0;  \
		gp ## id ##.e22 = 0.0;  \
		gp ## id ##.e12 = 0.0;  \
		gp ## id ##.E = E;      \
		gp ## id ##.niu = niu;  \
		gp ## id ##.Kf = Kf;    \
		gp ## id ##.k = k;      \
		gp ## id ##.miu = miu
		Init_Gauss_Point_Mat_Param(1);
		Init_Gauss_Point_Mat_Param(2);
		Init_Gauss_Point_Mat_Param(3);
		Init_Gauss_Point_Mat_Param(4);
	}
}

void Model_S2D_CHM_s_FEM_uUp::clear_mesh(void)
{
	if (nodes)
	{
		delete[] nodes;
		nodes = nullptr;
		node_x_num = 0;
		node_y_num = 0;
		node_num = 0;
	}
	if (elems)
	{
		delete[] elems;
		elems = nullptr;
		elem_x_num = 0;
		elem_y_num = 0;
		elem_num = 0;
	}
}
