#include "SimulationCore_pcp.h"

#include "Model_S2D_ME_s.h"

Model_S2D_ME_s::Model_S2D_ME_s() :
	nodes(nullptr), node_x_num(0), node_y_num(0), node_num(0),
	elems(nullptr), elem_x_num(0), elem_y_num(0), elem_num(0),
	pcls(nullptr), pcl_num(0),
	bfxs(nullptr), bfx_num(0), bfys(nullptr), bfy_num(0),
	txs(nullptr), tx_num(0), tys(nullptr), ty_num(0),
	axs(nullptr), ax_num(0), ays(nullptr), ay_num(0),
	vxs(nullptr), vx_num(0), vys(nullptr), vy_num(0) {}

Model_S2D_ME_s::~Model_S2D_ME_s()
{
	clear_mesh();
	clear_pcls();
	clear_bfxs();
	clear_bfys();
	clear_txs();
	clear_tys();
	clear_axs();
	clear_ays();
	clear_vxs();
	clear_vys();
}

int Model_S2D_ME_s::init_mesh(
	double _x0, double _y0, double _xn, double _yn,
	size_t _x_num, size_t _y_num, double h_tol_r)
{
	clear_mesh();
	if (_x0 == _xn || _y0 == _yn || _x_num == 0 || _y_num == 0)
		return -1;
	x0 = _x0;
	y0 = _y0;
	xn = _xn;
	yn = _yn;
	hx = (_xn - _x0) / double(_x_num);
	hy = (_yn - _y0) / double(_y_num);
	h_tol = h_tol_r * (hx + hy) * 0.5;
	// node
	node_x_num = _x_num + 1;
	node_y_num = _y_num + 1;
	node_num = node_x_num * node_y_num;
	nodes = new Node[node_num];
	size_t k;
	k = 0;
	for (size_t y_id = 0; y_id < node_y_num; ++y_id)
		for (size_t x_id = 0; x_id < node_x_num; ++x_id)
		{
			Node &n = nodes[k++];
			n.x_id = x_id;
			n.y_id = y_id;
		}
	// element
	elem_x_num = _x_num;
	elem_y_num = _y_num;
	elem_num = elem_x_num * elem_y_num;
	elems = new Element[elem_num];
	k = 0;
	for (size_t y_id = 0; y_id < elem_y_num; ++y_id)
		for (size_t x_id = 0; x_id < elem_x_num; ++x_id)
		{
			Element &e = elems[k++];
			e.x_id = x_id;
			e.y_id = y_id;
			e.pn1 = nodes + node_x_num * y_id + x_id;
			e.pn2 = e.pn1 + 1;
			e.pn4 = e.pn1 + node_x_num;
			e.pn3 = e.pn4 + 1;
		}
	elem_area = hx * hy;
	cal_shape_func_value(elem_sfv, 0.0, 0.0);
	return 0;
}

inline void Model_S2D_ME_s::clear_mesh()
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

void Model_S2D_ME_s::alloc_pcls(size_t num)
{
	clear_pcls();
	pcl_num = num;
	pcls = new Particle[num];
}

void Model_S2D_ME_s::init_pcls(
	size_t num,
	double m,
	double density,
	double _ar
	)
{
	alloc_pcls(num);
	for (size_t id = 0; id < num; ++id)
	{
		Particle &pcl = pcls[id];
		pcl.id = id;
		pcl.vx = 0.0;
		pcl.vy = 0.0;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
		pcl.m = m;
		pcl.density = density;
		pcl.ar = _ar;
		pcl.s11 = 0.0;
		pcl.s12 = 0.0;
		pcl.s22 = 0.0;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
	}
}

void Model_S2D_ME_s::init_pcls(
	double _x0, double _y0, double _xn, double _yn,
	size_t _x_num, size_t _y_num, double density)
{
	double pcl_hx = (_xn - _x0) / double(_x_num);
	double pcl_hy = (_yn - _y0) / double(_y_num);
	pcl_num = _x_num * _y_num;
	init_pcls(pcl_num, 0.0, density, pcl_hx / pcl_hy);
	double pcl_m = density * pcl_hx * pcl_hy;

	Particle *ppcl = pcls;
	for (size_t y_id = 0; y_id < _y_num; ++y_id)
	{
		for (size_t x_id = 0; x_id < _x_num; ++x_id)
		{
			ppcl[x_id].x = _x0 + pcl_hx * (0.5 + double(x_id));
			ppcl[x_id].y = _y0 + pcl_hy * (0.5 + double(y_id));
			ppcl[x_id].m = pcl_m;
		}
		ppcl += _x_num;
	}
}

inline void Model_S2D_ME_s::clear_pcls()
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
	}
	pcl_num = 0;
}

size_t Model_S2D_ME_s::distribute_pcl_mat_to_elems(Particle &pcl)
{
	// half length
	pcl.vol = pcl.m / pcl.density;
	double hlen = sqrt(pcl.vol / pcl.ar) * 0.5;
	double wlen = hlen * pcl.ar;
	double xl = pcl.x - wlen;
	double xu = pcl.x + wlen;
	double yl = pcl.y - hlen;
	double yu = pcl.y + hlen;

	// check if particles is out of mesh
	if (xu <= x0 || xl >= xn || yu <= y0 || yl >= yn)
		return 0;

	// trim particle if it lies at edge
	if (xl < x0)
		xl = x0;
	if (xu > xn)
		xu = xn;
	if (yl < y0)
		yl = y0;
	if (yu > yn)
		yu = yn;

	size_t xl_id = size_t(floor((xl - x0 + h_tol) / hx));
	size_t xu_id = size_t(ceil((xu - x0 - h_tol) / hx));
	size_t x_num = xu_id - xl_id;
	size_t yl_id = size_t(floor((yl - y0 + h_tol) / hy));
	size_t yu_id = size_t(ceil((yu - y0 - h_tol) / hy));
	size_t y_num = yu_id - yl_id;
	size_t elem_num = x_num * y_num;

	Element *pe;
	if (elem_num == 1)
	{
		pe = elems + elem_x_num * yl_id + xl_id;
		pe->ve_vol += (xu - xl) * (yu - yl);
		return 1;
	}

	double x_len1, x_len2, y_len1, y_len2;
	if (x_num == 1)
	{
		if (y_num == 2)
		{
			x_len1 = xu - xl;
			y_len1 = y0 + double(yl_id + 1) * hy - yl;
			y_len2 = yu - yl - y_len1;
			// elem 1
			pe = elems + elem_x_num * yl_id + xl_id;
			pe->ve_vol += x_len1 * y_len1;
			// elem 2
			pe += elem_x_num;
			pe->ve_vol += x_len1 * y_len2;
			return 2;
		}
	}
	else if (x_num == 2)
	{
		x_len1 = x0 + double(xl_id + 1) * hx - xl;
		x_len2 = xu - xl - x_len1;
		if (y_num == 1)
		{
			y_len1 = yu - yl;
			// elem 1
			pe = elems + elem_x_num * yl_id + xl_id;
			pe->ve_vol += x_len1 * y_len1;
			// elem 2
			++pe;
			pe->ve_vol += x_len2 * y_len1;
			return 2;
		}
		else if (y_num == 2)
		{
			y_len1 = y0 + double(yl_id + 1) * hy - yl;
			y_len2 = yu - yl - y_len1;
			// elem 1
			pe = elems + elem_x_num * yl_id + xl_id;
			pe->ve_vol += x_len1 * y_len1;
			// elem 2
			++pe;
			pe->ve_vol += x_len2 * y_len1;
			// elem 4
			pe += elem_x_num;
			pe->ve_vol += x_len2 * y_len2;
			// elem 3
			--pe;
			pe->ve_vol += x_len1 * y_len2;
			return 4;
		}
	}

	double *pcl_x_len, *pcl_y_len;
	// x
	pcl_x_len_buf.reset();
	if (x_num == 1)
	{
		pcl_x_len = &x_len1;
		x_len1 = xu - xl;
	}
	else
	{
		pcl_x_len = pcl_x_len_buf.alloc(x_num);
		pcl_x_len[0] = x0 + double(xl_id + 1) * hx - xl;
		for (size_t i = 1; i < x_num - 1; ++i)
			pcl_x_len[i] = hx;
		pcl_x_len[x_num - 1] = xu - double(xu_id - 1) * hx - x0;
	}
	// y
	pcl_y_len_buf.reset();
	if (y_num == 1)
	{
		pcl_y_len = &y_len1;
		y_len1 = yu - yl;
	}
	else
	{
		pcl_y_len = pcl_y_len_buf.alloc(y_num);
		pcl_y_len[0] = y0 + double(yl_id + 1) * hy - yl;
		for (size_t i = 1; i < y_num - 1; ++i)
			pcl_y_len[i] = hy;
		pcl_y_len[y_num - 1] = yu - double(yu_id - 1) * hy - y0;
	}

	pe = elems + elem_x_num * yl_id + xl_id;
	for (size_t y_id = 0; y_id < y_num; ++y_id)
	{
		for (size_t x_id = 0; x_id < x_num; ++x_id)
			pe[x_id].ve_vol += pcl_x_len[x_id] * pcl_y_len[y_id];
		pe += elem_x_num;
	}

	return elem_num;
}
