#include "SimulationCore_pcp.h"

#include "Model_S2D_ME_s_up.h"

Model_S2D_ME_s_up::Model_S2D_ME_s_up() :
	Model("Model_S2D_ME_s_up"),
	elems(nullptr), elem_x_num(0), elem_y_num(0), elem_num(0),
	nodes(nullptr), node_x_num(0), node_y_num(0), node_num(0),
	pcls(nullptr), pcl_num(0),
	bfxs(nullptr), bfx_num(0),
	bfys(nullptr), bfy_num(0),
	txs(nullptr), tx_num(0),
	tys(nullptr), ty_num(0),
	uxs(nullptr), ux_num(0),
	uys(nullptr), uy_num(0), 
	pbcs(nullptr), pbc_num(0) {}

Model_S2D_ME_s_up::~Model_S2D_ME_s_up()
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
	if (pbcs)
	{
		delete[] pbcs;
		pbcs = nullptr;
		pbc_num = 0;
	}
}

void Model_S2D_ME_s_up::init_pcl(
	size_t num, double m, double density, double _E, double _niu, double _K)
{
	clear_pcl();
	pcl_num = num;
	pcls = new Particle[pcl_num];
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.index = pcl_id;
		pcl.m = m;
		pcl.density = density;
		pcl.ax = 0.0;
		pcl.ay = 0.0;
		pcl.vx = 0.0;
		pcl.vy = 0.0;
		pcl.s11 = 0.0;
		pcl.s22 = 0.0;
		pcl.s12 = 0.0;
		pcl.p = 0.0;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
	}
	E = _E;
	niu = _niu;
	K = _K;
}

void Model_S2D_ME_s_up::clear_pcl(void)
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
		pcl_num = 0;
	}
}

void Model_S2D_ME_s_up::init_mesh(
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

	elem_vol = h * h;
	dof_num = node_num * 3;
	gp_w = elem_vol * 0.25;
	cal_shape_func_value(sf1, -0.5773502692, -0.5773502692);
	cal_shape_func_value(sf2, 0.5773502692, -0.5773502692);
	cal_shape_func_value(sf3, 0.5773502692, 0.5773502692);
	cal_shape_func_value(sf4, -0.5773502692, 0.5773502692);
#define Form_dN_dx_mat(id)                              \
	dN_dx_mat ## id ## [0][0] = sf ## id ## .dN1_dx; \
	dN_dx_mat ## id ## [0][1] = sf ## id ## .dN2_dx; \
	dN_dx_mat ## id ## [0][2] = sf ## id ## .dN3_dx; \
	dN_dx_mat ## id ## [0][3] = sf ## id ## .dN4_dx; \
	dN_dx_mat ## id ## [0][4] = 0.0;                    \
	dN_dx_mat ## id ## [0][5] = 0.0;                    \
	dN_dx_mat ## id ## [0][6] = 0.0;                    \
	dN_dx_mat ## id ## [0][7] = 0.0;                    \
	dN_dx_mat ## id ## [1][0] = 0.0;                    \
	dN_dx_mat ## id ## [1][1] = 0.0;                    \
	dN_dx_mat ## id ## [1][2] = 0.0;                    \
	dN_dx_mat ## id ## [1][3] = 0.0;                    \
	dN_dx_mat ## id ## [1][4] = sf ## id ## .dN1_dy; \
	dN_dx_mat ## id ## [1][5] = sf ## id ## .dN2_dy; \
	dN_dx_mat ## id ## [1][6] = sf ## id ##.dN3_dy; \
	dN_dx_mat ## id ## [1][7] = sf ## id ##.dN4_dy; \
	dN_dx_mat ## id ## [2][0] = sf ## id ##.dN1_dy; \
	dN_dx_mat ## id ## [2][1] = sf ## id ##.dN2_dy; \
	dN_dx_mat ## id ## [2][2] = sf ## id ##.dN3_dy; \
	dN_dx_mat ## id ## [2][3] = sf ## id ##.dN4_dy; \
	dN_dx_mat ## id ## [2][4] = sf ## id ##.dN1_dx; \
	dN_dx_mat ## id ## [2][5] = sf ## id ##.dN2_dx; \
	dN_dx_mat ## id ## [2][6] = sf ## id ##.dN3_dx; \
	dN_dx_mat ## id ## [2][7] = sf ## id ##.dN4_dx
	Form_dN_dx_mat(1);
	Form_dN_dx_mat(2);
	Form_dN_dx_mat(3);
	Form_dN_dx_mat(4);
}

void Model_S2D_ME_s_up::clear_mesh(void)
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
