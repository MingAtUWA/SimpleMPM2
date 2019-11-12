#include "SimulationCore_pcp.h"

#include "Model_S2D_CHM_s_uUp.h"

Model_S2D_CHM_s_uUp::Model_S2D_CHM_s_uUp() :
	Model("Model_S2D_CHM_s_uUp"),
	nodes(nullptr), node_x_num(0), node_y_num(0), node_num(0),
	elems(nullptr), elem_x_num(0), elem_y_num(0), elem_num(0),
	pcls(nullptr), pcl_num(0),
	bfxs(nullptr), bfx_num(0), bfys(nullptr), bfy_num(0),
	txs(nullptr),  tx_num(0),  tys(nullptr),  ty_num(0),
	usxs(nullptr), usx_num(0), usys(nullptr), usy_num(0),
	ufxs(nullptr), ufx_num(0), ufys(nullptr), ufy_num(0),
	pbcs(nullptr), pbc_num(0) {}

Model_S2D_CHM_s_uUp::~Model_S2D_CHM_s_uUp()
{
	clear_mesh();
	clear_pcl();
	if (bfxs)
	{
		delete[] bfxs;
		bfxs = nullptr;
	}
	bfx_num = 0;
	if (bfys)
	{
		delete[] bfys;
		bfys = nullptr;
	}
	bfy_num = 0;
	if (txs)
	{
		delete[] txs;
		txs = nullptr;
	}
	tx_num = 0;
	if (tys)
	{
		delete[] tys;
		tys = nullptr;
	}
	ty_num = 0;
	// solid bc
	if (usxs)
	{
		delete[] usxs;
		usxs = nullptr;
	}
	usx_num = 0;
	if (usys)
	{
		delete[] usys;
		usys = nullptr;
	}
	usy_num = 0;
	// fluid bc
	if (ufxs)
	{
		delete[] ufxs;
		ufxs = nullptr;
	}
	ufx_num = 0;
	if (ufys)
	{
		delete[] ufys;
		ufys = nullptr;
	}
	ufy_num = 0;
}

void Model_S2D_CHM_s_uUp::init_mesh(
	double _h, size_t _elem_x_num, size_t _elem_y_num,
	double x_start, double y_start, double h_tol_ratio)
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
}

void Model_S2D_CHM_s_uUp::clear_mesh(void)
{
	if (nodes)
	{
		delete[] nodes;
		nodes = nullptr;
	}
	node_x_num = 0;
	node_y_num = 0;
	node_num = 0;
	if (elems)
	{
		delete[] elems;
		elems = nullptr;
	}
	elem_x_num = 0;
	elem_y_num = 0;
	elem_num = 0;
}

void Model_S2D_CHM_s_uUp::init_pcl(size_t num,
	double n, double m_s, double density_s, double density_f,
	double E, double niu, double Kf, double k, double miu)
{
	clear_pcl();
	pcl_num = num;
	pcls = new Particle[num];
	for (size_t pcl_id = 0; pcl_id < num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.index = pcl_id;
		pcl.n = n;
		pcl.m_s = m_s;
		pcl.density_s = density_s;
		pcl.density_f = density_f;
		pcl.ux_s = 0.0;
		pcl.uy_s = 0.0;
		pcl.vx_s = 0.0;
		pcl.vy_s = 0.0;
		pcl.ax_s = 0.0;
		pcl.ay_s = 0.0;
		pcl.ux_f = 0.0;
		pcl.uy_f = 0.0;
		pcl.vx_f = 0.0;
		pcl.vy_f = 0.0;
		pcl.ax_f = 0.0;
		pcl.ay_f = 0.0;
		pcl.s11 = 0.0;
		pcl.s12 = 0.0;
		pcl.s22 = 0.0;
		pcl.p = 0.0;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
		pcl.E = E;
		pcl.niu = niu;
		pcl.Kf = Kf;
		pcl.k = k;
		pcl.miu = miu;
	}
}

void Model_S2D_CHM_s_uUp::clear_pcl(void)
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
	}
	pcl_num = 0;
}
