#include "SimulationCore_pcp.h"

#include "Model_S2D_ME_s_FEM_up.h"

Model_S2D_ME_s_FEM_up::Model_S2D_ME_s_FEM_up() :
	Model("Model_S2D_ME_s_FEM_up"),
	elems(nullptr), elem_x_num(0), elem_y_num(0), elem_num(0),
	nodes(nullptr), node_x_num(0), node_y_num(0), node_num(0),
	//bfxs(nullptr), bfx_num(0),
	//bfys(nullptr), bfy_num(0),
	txs(nullptr), tx_num(0),
	tys(nullptr), ty_num(0),
	//pbcs(nullptr), pbc_num(0),
	uxs(nullptr), ux_num(0),
	uys(nullptr), uy_num(0) {}

Model_S2D_ME_s_FEM_up::~Model_S2D_ME_s_FEM_up()
{
	clear_mesh();
	//if (bfxs)
	//{
	//	delete[] bfxs;
	//	bfxs = nullptr;
	//	bfx_num = 0;
	//}
	//if (bfys)
	//{
	//	delete[] bfys;
	//	bfys = nullptr;
	//	bfy_num = 0;
	//}
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
	//if (pbcs)
	//{
	//	delete[] pbcs;
	//	pbcs = nullptr;
	//	pbc_num = 0;
	//}
	if (uxs)
	{
		delete[] uxs;
		uxs = nullptr;
		ux_num = 0;
	}
	if (uys)
	{
		delete[] uys;
		uys = nullptr;
		uy_num = 0;
	}
}

void Model_S2D_ME_s_FEM_up::init_mesh(
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

	dof_num = node_num * 3;
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

void Model_S2D_ME_s_FEM_up::clear_mesh(void)
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

void Model_S2D_ME_s_FEM_up::init_mat_param(
	double density, double E, double niu, double K)
{
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		n.ux = 0.0;
		n.uy = 0.0;
		n.p = 0.0;
	}
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &e = elems[e_id];
#define Init_Gauss_Point_Mat_Param(id)        \
		GaussPoint &gp ## id = e.gp ## id ##; \
		gp ## id ##.density = density;        \
		gp ## id ##.E = E;     \
		gp ## id ##.niu = niu; \
		gp ## id ##.K = K;     \
		gp ## id ##.ax = 0.0;  \
		gp ## id ##.ay = 0.0;  \
		gp ## id ##.vx = 0.0;  \
		gp ## id ##.vy = 0.0;  \
		gp ## id ##.ux = 0.0;  \
		gp ## id ##.uy = 0.0;  \
		gp ## id ##.s11 = 0.0; \
		gp ## id ##.s22 = 0.0; \
		gp ## id ##.s12 = 0.0; \
		gp ## id ##.p = 0.0;   \
		gp ## id ##.e11 = 0.0; \
		gp ## id ##.e22 = 0.0; \
		gp ## id ##.e12 = 0.0
		Init_Gauss_Point_Mat_Param(1);
		Init_Gauss_Point_Mat_Param(2);
		Init_Gauss_Point_Mat_Param(3);
		Init_Gauss_Point_Mat_Param(4);
	}
}
