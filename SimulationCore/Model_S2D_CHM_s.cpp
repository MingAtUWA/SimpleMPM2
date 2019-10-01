#include "SimulationCore_pcp.h"

#include "Model_S2D_CHM_s.h"

Model_S2D_CHM_s::Model_S2D_CHM_s() :
	Model("Model_S2D_CHM_s"),
	elems(nullptr), elem_x_num(0), elem_y_num(0), elem_num(0),
	nodes(nullptr), node_x_num(0), node_y_num(0), node_num(0),
	pcls(nullptr), pcl_num(0),
	bfx_num(0), bfy_num(0), bfxs(nullptr), bfys(nullptr),
	tx_num(0),  ty_num(0),  txs(nullptr),  tys(nullptr),
	asx_num(0), asy_num(0), asxs(nullptr), asys(nullptr),
	vsx_num(0), vsy_num(0), vsxs(nullptr), vsys(nullptr),
	afx_num(0), afy_num(0), afxs(nullptr), afys(nullptr),
	vfx_num(0), vfy_num(0), vfxs(nullptr), vfys(nullptr),
	pcl_var_mem(50), x_var_info_buf(20), y_var_info_buf(20) {}

Model_S2D_CHM_s::~Model_S2D_CHM_s()
{
	clear_mesh();
	clear_pcl();
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
	// solid bc
	if (asxs)
	{
		delete[] asxs;
		asxs = nullptr;
		asx_num = 0;
	}
	if (asys)
	{
		delete[] asys;
		asys = nullptr;
		asy_num = 0;
	}
	if (vsxs)
	{
		delete[] vsxs;
		vsxs = nullptr;
		vsx_num = 0;
	}
	if (vsys)
	{
		delete[] vsys;
		vsys = nullptr;
		vsy_num = 0;
	}
	// fluid bc
	if (afxs)
	{
		delete[] afxs;
		afxs = nullptr;
		afx_num = 0;
	}
	if (afys)
	{
		delete[] afys;
		afys = nullptr;
		afy_num = 0;
	}
	if (vfxs)
	{
		delete[] vfxs;
		vfxs = nullptr;
		vfx_num = 0;
	}
	if (vfys)
	{
		delete[] vfys;
		vfys = nullptr;
		vfy_num = 0;
	}
	pcl_var_mem.clear();
	x_var_info_buf.clear();
	y_var_info_buf.clear();
}

void Model_S2D_CHM_s::init_mesh(
	double _h, size_t _elem_x_num, size_t _elem_y_num,
	double x_start, double y_start)
{
	clear_mesh();
	h = _h;
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
			Element &elem = elems[k];
			elem.index_x = j;
			elem.index_y = i;
			++k;
		}
	k = 0;
	for (i = 0; i < node_y_num; ++i)
		for (j = 0; j < node_x_num; ++j)
		{
			Node &node = nodes[k];
			node.index_x = j;
			node.index_y = i;
			++k;
		}
}

void Model_S2D_CHM_s::clear_mesh(void)
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

void Model_S2D_CHM_s::init_pcl(size_t num,
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
		pcl.vx_s = 0.0;
		pcl.vy_s = 0.0;
		pcl.vx_f = 0.0;
		pcl.vy_f = 0.0;
		pcl.ux_s = 0.0;
		pcl.uy_s = 0.0;
		pcl.ux_f = 0.0;
		pcl.uy_f = 0.0;
		pcl.n = n;
		pcl.m_s = m_s;
		pcl.density_s = density_s;
		pcl.density_f = density_f;
		pcl.E = E;
		pcl.niu = niu;
		pcl.Kf = Kf;
		pcl.k = k;
		pcl.miu = miu;
		pcl.s11 = 0.0;
		pcl.s12 = 0.0;
		pcl.s22 = 0.0;
		pcl.p = 0.0;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
	}
	pcl_var_mem.reset();
	pcl_var_mem.set_page_size(pcl_num);
}

void Model_S2D_CHM_s::clear_pcl(void)
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
		pcl_num = 0;
	}
	pcl_var_mem.clear();
}


bool Model_S2D_CHM_s::init_pcl_standard(Particle &pcl)
{
	// check if particles is out of mesh
	if (pcl.x < x0 || pcl.x >= xn || pcl.y < y0 || pcl.y >= yn)
	{
		pcl.elem_num = 0;
		return false;
	}

	ParticleCalVar &pcl_var = pcl.var;
	pcl_var.vol = pcl.m_s / (pcl.density_s * (1.0 - pcl.n));
	pcl_var.elem_x_id = size_t((pcl.x - x0) / h);
	pcl_var.elem_y_id = size_t((pcl.y - y0) / h);
	cal_pcl_shape_func_value(pcl_var);

	return true;
}

bool Model_S2D_CHM_s::init_pcl_GIMP(Particle &pcl)
{
	ParticleCalVar &pcl_var = pcl.var;
	pcl_var.vol = pcl.m_s / ((1.0 - pcl.n) * pcl.density_s);

	double hlen = sqrt(pcl_var.vol) * 0.5; // half length
	double xl = pcl.x - hlen;
	double xu = pcl.x + hlen;
	double yl = pcl.y - hlen;
	double yu = pcl.y + hlen;

	// check if particles is out of mesh
	if (xu < x0 || xl >= xn || yu < y0 || yl >= yn)
	{
		pcl.elem_num = 0;
		pcl.vars = nullptr;
		return false;
	}

	// == 00000001b at external edge
	// == 00000010b at internal edge
	// == 00000011b at both external and internal edge
	unsigned char pos_tag = 0;

	// trim particle if it lies at the edge
	if (xl < x0)
	{
		xl = x0;
		pos_tag = 1;
	}
	if (xu > xn)
	{
		xu = xn;
		pos_tag = 1;
	}
	if (yl < y0)
	{
		yl = y0;
		pos_tag = 1;
	}
	if (yu > yn)
	{
		yu = yn;
		pos_tag = 1;
	}

	size_t xl_id = size_t((xl - x0) / h);
	size_t xu_id = size_t((xu - x0) / h);
	if (xu - x0 > h * double(xu_id)) ++xu_id;
	size_t x_num = xu_id - xl_id;
	size_t yl_id = size_t((yl - y0) / h);
	size_t yu_id = size_t((yu - y0) / h);
	if (yu - y0 > h * double(yu_id)) ++yu_id;
	size_t y_num = yu_id - yl_id;
	pcl.elem_num = x_num * y_num;

	if (pos_tag)
	{
		// part of this particle is outside mesh
		pcl_var.x = (xl + xu) * 0.5;
		pcl_var.y = (yl + yu) * 0.5;
		pcl_var.vol = (xu - xl) * (yu - yl);
	}
	else
	{
		// all particle is inside mesh
		pcl_var.x = pcl.x;
		pcl_var.y = pcl.y;
	}
	pcl_var.elem_x_id = xl_id;
	pcl_var.elem_y_id = yl_id;
	cal_pcl_shape_func_value(pcl_var);

	if (pcl.elem_num == 1)
	{
		pcl.vars = &pcl.var;
		return true;
	}

	// particle is at internal edge
	//pos_tag |= 2;
	pcl.vars = pcl_var_mem.alloc_fast(pcl.elem_num);
	double x_len1, x_len2, y_len1, y_len2;
	if (x_num == 1)
	{
		if (y_num == 2)
		{
			x_len1 = xu - xl;
			y_len1 = y0 + double(yl_id + 1) * h - yl;
			y_len2 = yu - yl - y_len1;
			ParticleCalVar &pcl_var1 = pcl.vars[0];
			pcl_var1.x = pcl.x;
			pcl_var1.y = yl + y_len1 * 0.5;
			pcl_var1.vol = x_len1 * y_len1;
			pcl_var1.elem_x_id = xl_id;
			pcl_var1.elem_y_id = yl_id;
			cal_pcl_shape_func_value(pcl_var1);
			ParticleCalVar &pcl_var2 = pcl.vars[1];
			pcl_var2.x = pcl.x;
			pcl_var2.y = yu - y_len2 * 0.5;
			pcl_var2.vol = x_len1 * y_len2;
			pcl_var2.elem_x_id = xl_id;
			pcl_var2.elem_y_id = yl_id + 1;
			cal_pcl_shape_func_value(pcl_var2);
			return true;
		}
	}
	else if (x_num == 2)
	{
		x_len1 = x0 + double(xl_id + 1) * h - xl;
		x_len2 = xu - xl - x_len1;
		if (y_num == 1)
		{
			y_len1 = yu - yl;
			ParticleCalVar &pcl_var1 = pcl.vars[0];
			pcl_var1.x = xl + x_len1 * 0.5;
			pcl_var1.y = pcl.y;
			pcl_var1.vol = x_len1 * y_len1;
			pcl_var1.elem_x_id = xl_id;
			pcl_var1.elem_y_id = yl_id;
			cal_pcl_shape_func_value(pcl_var1);
			ParticleCalVar &pcl_var2 = pcl.vars[1];
			pcl_var2.x = xu - x_len2 * 0.5;
			pcl_var2.y = pcl.y;
			pcl_var2.vol = x_len2 * y_len1;
			pcl_var2.elem_x_id = xl_id + 1;
			pcl_var2.elem_y_id = yl_id;
			cal_pcl_shape_func_value(pcl_var2);
			return true;
		}
		else
		{
			y_len1 = y0 + double(yl_id + 1) * h - yl;
			y_len2 = yu - yl - y_len1;
			ParticleCalVar &pcl_var1 = pcl.vars[0];
			pcl_var1.x = xl + x_len1 * 0.5;
			pcl_var1.y = yl + y_len1 * 0.5;
			pcl_var1.vol = x_len1 * y_len1;
			pcl_var1.elem_x_id = xl_id;
			pcl_var1.elem_y_id = yl_id;
			cal_pcl_shape_func_value(pcl_var1);
			ParticleCalVar &pcl_var2 = pcl.vars[1];
			pcl_var2.x = xu - x_len2 * 0.5;
			pcl_var2.y = yl + y_len1 * 0.5;
			pcl_var2.vol = x_len2 * y_len1;
			pcl_var2.elem_x_id = xl_id + 1;
			pcl_var2.elem_y_id = yl_id;
			cal_pcl_shape_func_value(pcl_var2);
			ParticleCalVar &pcl_var3 = pcl.vars[2];
			pcl_var3.x = xl + x_len1 * 0.5;
			pcl_var3.y = yu - y_len2 * 0.5;
			pcl_var3.vol = x_len1 * y_len2;
			pcl_var3.elem_x_id = xl_id;
			pcl_var3.elem_y_id = yl_id + 1;
			cal_pcl_shape_func_value(pcl_var3);
			ParticleCalVar &pcl_var4 = pcl.vars[3];
			pcl_var4.x = xu - x_len2 * 0.5;
			pcl_var4.y = yu - y_len2 * 0.5;
			pcl_var4.vol = x_len2 * y_len2;
			pcl_var4.elem_x_id = xl_id + 1;
			pcl_var4.elem_y_id = yl_id + 1;
			cal_pcl_shape_func_value(pcl_var4);
			return true;
		}
	}

	PclVarInfo *x_var_infos, *y_var_infos;
	// x
	x_var_info_buf.reset();
	if (x_num == 1)
	{
		x_var_infos = x_var_info_buf.alloc();
		x_var_infos[0].len = xu - xl;
		x_var_infos[0].pos = (xu + xl) * 0.5;
		x_var_infos[0].elem_id = size_t((x_var_infos[0].pos - x0) / h);
	}
	else
	{
		x_var_infos = x_var_info_buf.alloc(x_num);
		x_var_infos[0].len = x0 + double(xl_id + 1) * h - xl;
		x_var_infos[0].pos = xl + x_var_infos[0].len * 0.5;
		x_var_infos[0].elem_id = size_t((x_var_infos[0].pos - x0) / h);
		for (size_t i = 1; i < x_num - 1; ++i)
		{
			PclVarInfo &var_info = x_var_infos[i];
			var_info.len = h;
			var_info.pos = x_var_infos[i - 1].pos + (x_var_infos[i - 1].len + h) * 0.5;
			var_info.elem_id = x_var_infos[i - 1].elem_id + 1;
		}
		x_var_infos[x_num - 1].len = xu - double(xu_id - 1) * h - x0;
		x_var_infos[x_num - 1].pos = xu - x_var_infos[x_num - 1].len * 0.5;
		x_var_infos[x_num - 1].elem_id = x_var_infos[x_num - 2].elem_id + 1;
	}
	// y
	y_var_info_buf.reset();
	if (y_num == 1)
	{
		y_var_infos = y_var_info_buf.alloc();
		y_var_infos[0].len = yu - yl;
		y_var_infos[0].pos = (yu + yl) * 0.5;
		y_var_infos[0].elem_id = size_t((y_var_infos[0].pos - y0) / h);
	}
	else
	{
		y_var_infos = y_var_info_buf.alloc(y_num);
		y_var_infos[0].len = y0 + double(yl_id + 1) * h - yl;
		y_var_infos[0].pos = yl + y_var_infos[0].len * 0.5;
		y_var_infos[0].elem_id = size_t((y_var_infos[0].pos - y0) / h);
		for (size_t i = 1; i < y_num - 1; ++i)
		{
			PclVarInfo &var_info = y_var_infos[i];
			var_info.len = h;
			var_info.pos = y_var_infos[i - 1].pos + (y_var_infos[i - 1].len + h) * 0.5;
			var_info.elem_id = y_var_infos[i - 1].elem_id + 1;
		}
		y_var_infos[y_num - 1].len = yu - double(yu_id - 1) * h - x0;
		y_var_infos[y_num - 1].pos = yu - y_var_infos[y_num - 1].len * 0.5;
		y_var_infos[y_num - 1].elem_id = size_t((y_var_infos[y_num - 1].pos - y0) / h);
	}

	size_t k = 0;
	for (size_t j = 0; j < y_num; ++j)
		for (size_t i = 0; i < x_num; ++i)
		{
			ParticleCalVar &pcl_var = pcl.vars[k];
			pcl_var.x = x_var_infos[i].pos;
			pcl_var.y = y_var_infos[j].pos;
			pcl_var.vol = x_var_infos[i].len * y_var_infos[j].len;
			pcl_var.elem_x_id = x_var_infos[i].elem_id;
			pcl_var.elem_y_id = y_var_infos[j].elem_id;
			cal_pcl_shape_func_value(pcl_var);
			++k;
		}

	return true;
}
