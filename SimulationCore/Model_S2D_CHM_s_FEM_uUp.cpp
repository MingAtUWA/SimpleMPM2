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
	elem_vol = h * h;
	x0 = x_start;
	xn = x0 + _h * double(_elem_x_num);
	y0 = y_start;
	yn = y0 + _h * double(_elem_y_num);
	elem_x_num = _elem_x_num;
	elem_y_num = _elem_y_num;
	elem_num = elem_x_num * elem_y_num;
	elems = new Element[elem_num];
	node_x_num = _elem_x_num + 1;
	node_y_num = _elem_y_num + 1;
	node_num = node_x_num * node_y_num;
	nodes = new Node[node_num];
	size_t i, j, k;
	k = 0;
	for (i = 0; i < elem_y_num; ++i)
		for (j = 0; j < elem_x_num; ++j)
		{
			Element &e = elems[k++];
			e.index_x = j;
			e.index_y = i;
			e.e11 = 0.0;
			e.e22 = 0.0;
			e.e12 = 0.0;
			e.s11 = 0.0;
			e.s22 = 0.0;
			e.s12 = 0.0;
			e.p = 0.0;
		}
	k = 0;
	for (i = 0; i < node_y_num; ++i)
		for (j = 0; j < node_x_num; ++j)
		{
			Node &n = nodes[k++];
			n.index_x = j;
			n.index_y = i;
			n.x = x0 + h * double(j);
			n.y = y0 + h * double(i);
			n.p = 0.0;
			// solid phase
			n.ax_s = 0.0;
			n.ay_s = 0.0;
			n.vx_s = 0.0;
			n.vy_s = 0.0;
			n.ux_s = 0.0;
			n.uy_s = 0.0;
			// for fluid phase
			n.ax_f = 0.0;
			n.ay_f = 0.0;
			n.vx_f = 0.0;
			n.vy_f = 0.0;
			n.ux_f = 0.0;
			n.uy_f = 0.0;
		}
}

void Model_S2D_CHM_s_FEM_uUp::init_mat_param(
	double n, double density_s, double density_f,
	double E, double niu, double Kf, double k, double miu)
{
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &e = elems[e_id];
		e.n = n;
		e.density_s = density_s;
		e.density_f = density_f;
		e.E = E;
		e.niu = niu;
		e.Kf = Kf;
		e.k = k;
		e.miu = miu;
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
