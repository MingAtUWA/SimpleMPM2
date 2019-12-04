#include "SimulationCore_pcp.h"

#include "Model_T2D_CHM_s.h"

Model_T2D_CHM_s::Model_T2D_CHM_s() :
	Model("Model_T2D_CHM_s"),
	elems(nullptr), elem_num(0), nodes(nullptr), node_num(0),
	pcls(nullptr), pcl_num(0),
	bfx_num(0), bfy_num(0), bfxs(nullptr), bfys(nullptr),
	tx_num(0),  ty_num(0),  txs(nullptr),  tys(nullptr),
	asx_num(0), asy_num(0), asxs(nullptr), asys(nullptr),
	vsx_num(0), vsy_num(0), vsxs(nullptr), vsys(nullptr),
	afx_num(0), afy_num(0), afxs(nullptr), afys(nullptr),
	vfx_num(0), vfy_num(0), vfxs(nullptr), vfys(nullptr),
	grid_x_num(0), grid_y_num(0), grid_num(0), bg_grids(nullptr) {}

Model_T2D_CHM_s::~Model_T2D_CHM_s()
{
	clear_mesh();
	clear_pcls();
	clear_bfxs();
	clear_bfys();
	clear_txs();
	clear_tys();
	clear_asxs();
	clear_asys();
	clear_afxs();
	clear_afys();
	clear_vsxs();
	clear_vsys();
	clear_vfxs();
	clear_vfys();
	clear_bg_mesh();
}

void Model_T2D_CHM_s::init_mesh(double *node_coords, size_t n_num,
								size_t *elem_n_ids,  size_t e_num)
{
	clear_mesh();
	if (n_num == 0 || e_num == 0)
		return;
	
	node_num = n_num;
	nodes = new Node[node_num];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		n.id = n_id;
		n.x = node_coords[2*n_id];
		n.y = node_coords[2*n_id + 1];
	}
	
	elem_num = e_num;
	elems = new Element[elem_num];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &e = elems[e_id];
		e.id = e_id;
		e.n1 = elem_n_ids[3*e_id];
		e.n2 = elem_n_ids[3*e_id + 1];
		e.n3 = elem_n_ids[3*e_id + 2];
		Node &n1 = nodes[e.n1];
		Node &n2 = nodes[e.n2];
		Node &n3 = nodes[e.n3];
		e.area_2 = (n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (n2.y - n1.y);
		// Ensure that nodes index of element is in counter-clockwise sequence.
		if (e.area_2 < 0.0)
		{
			e.area_2 = -e.area_2;
			size_t n_tmp = e.n2;
			e.n2 = e.n3;
			e.n3 = n_tmp;
		}
		e.area = e.area_2 * 0.5;
		Node &_n2 = nodes[e.n2];
		Node &_n3 = nodes[e.n3];
		ShapeFuncValue &sf = e.sf;
		sf.dN1_dx = (_n2.y - _n3.y) / e.area_2;
		sf.dN1_dy = (_n3.x - _n2.x) / e.area_2;
		sf.dN2_dx = (_n3.y -  n1.y) / e.area_2;
		sf.dN2_dy = ( n1.x - _n3.x) / e.area_2;
		sf.dN3_dx = ( n1.y - _n2.y) / e.area_2;
		sf.dN3_dy = (_n2.x -  n1.x) / e.area_2;
		//e.sf.cal_shape_func_value(1.0/3.0, 1.0/3.0, n1.x, n1.y,
		//						  _n2.x, _n2.y, _n3.x, _n3.y);
	}
}

void Model_T2D_CHM_s::clear_mesh(void)
{
	if (nodes)
	{
		delete[] nodes;
		nodes = nullptr;
	}
	node_num = 0;
	if (elems)
	{
		delete[] elems;
		elems = nullptr;
	}
	elem_num = 0;
}

void Model_T2D_CHM_s::init_pcls(size_t num,
	double n, double m_s, double density_s, double density_f,
	double _E, double _niu, double _Kf, double _k, double _miu)
{
	clear_pcls();
	pcl_num = num;
	pcls = new Particle[pcl_num];
	for (size_t pcl_id = 0; pcl_id < num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.id = pcl_id;
		pcl.ux_s = 0.0;
		pcl.uy_s = 0.0;
		pcl.vx_s = 0.0;
		pcl.vy_s = 0.0;
		pcl.vx_f = 0.0;
		pcl.vy_f = 0.0;
		pcl.ux_f = 0.0;
		pcl.uy_f = 0.0;
		pcl.n = n;
		pcl.m_s = m_s;
		pcl.density_s = density_s;
		pcl.density_f = density_f;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
		pcl.s11 = 0.0;
		pcl.s12 = 0.0;
		pcl.s22 = 0.0;
		pcl.p = 0.0;
	}
	// constitutive model
	E = _E;
	niu = _niu;
	Kf = _Kf;
	k = _k;
	miu = _miu;
}

void Model_T2D_CHM_s::clear_pcls(void)
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
	}
	pcl_num = 0;
}

void Model_T2D_CHM_s::init_mesh(TriangleMesh &tri_mesh)
{
	clear_mesh();
	node_num = tri_mesh.get_node_num();
	elem_num = tri_mesh.get_elem_num();
	if (node_num == 0 || elem_num == 0)
		return;

	TriangleMesh::Node *ext_nodes = tri_mesh.get_nodes();
	nodes = new Node[node_num];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		TriangleMesh::Node &ext_n = ext_nodes[n_id];
		Node &n = nodes[n_id];
		n.id = ext_n.id;
		n.x  = ext_n.x;
		n.y  = ext_n.y;
	}

	TriangleMesh::Element *ext_elems = tri_mesh.get_elems();
	elems = new Element[elem_num];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		TriangleMesh::Element &ext_e = ext_elems[e_id];
		Element &e = elems[e_id];
		e.id = ext_e.id;
		e.n1 = ext_e.n1;
		e.n2 = ext_e.n2;
		e.n3 = ext_e.n3;
		e.area_2 = ext_e.area;
		e.area = e.area_2 * 0.5;
		Node &n1 = nodes[e.n1];
		Node &n2 = nodes[e.n2];
		Node &n3 = nodes[e.n3];
		ShapeFuncValue &sf = e.sf;
		sf.dN1_dx = (n2.y - n3.y) / e.area_2;
		sf.dN1_dy = (n3.x - n2.x) / e.area_2;
		sf.dN2_dx = (n3.y - n1.y) / e.area_2;
		sf.dN2_dy = (n1.x - n3.x) / e.area_2;
		sf.dN3_dx = (n1.y - n2.y) / e.area_2;
		sf.dN3_dy = (n2.x - n1.x) / e.area_2;
	}
}

void Model_T2D_CHM_s::init_pcls(TriangleMeshToParticles &mh_2_pcl,
	double n, double density_s, double density_f,
	double E, double niu, double Kf, double k, double miu)
{
	pcl_num = mh_2_pcl.get_pcl_num();
	if (pcl_num == 0)
		return;
	init_pcls(pcl_num, n, (1.0-n)*density_s, density_s, density_f, E, niu, Kf, k, miu);
	TriangleMeshToParticles::Particle *ext_ppcl = mh_2_pcl.first();
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.x = ext_ppcl->x;
		pcl.y = ext_ppcl->y;
		pcl.m_s *= ext_ppcl->vol;
		ext_ppcl = mh_2_pcl.next(ext_ppcl);
	}
}

int Model_T2D_CHM_s::apply_rigid_body_to_bg_mesh(double dtime)
{
	// reset reaction force
	rigid_circle.rfx = 0.0;
	rigid_circle.rfy = 0.0;
	rigid_circle.rm = 0.0;
	// update state
	rigid_circle.update(dtime);

	for (size_t pcl_id = 0; pcl_id < rigid_circle.pcl_num; ++pcl_id)
	{
		RCParticle &pcl = rigid_circle.pcls[pcl_id];
		apply_pcl_to_mesh(pcl);
	}
	
	double fx, fy;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		if (n.has_rb && n.has_mp && n.vol_rb != 0.0)
		{
			n.vx_rb /= n.vol_rb;
			n.vy_rb /= n.vol_rb;
			// rough contact for both solid and fluid phase
			fx = -((n.vx_rb - n.vx_s) * n.m_s + (n.vx_rb - n.vx_f) * n.m_f) / dtime;
			fy = -((n.vy_rb - n.vy_s) * n.m_s + (n.vy_rb - n.vy_f) * n.m_f) / dtime;
			rigid_circle.add_reaction_force(n.x, n.y, fx, fy);
			// modify nodal velocity
			n.vx_s = n.vx_rb;
			n.vy_s = n.vy_rb;
			n.vx_f = n.vx_rb;
			n.vy_f = n.vy_rb;
		}
	}

	return 0;
}

int Model_T2D_CHM_s::init_bg_mesh(double hx, double hy)
{
	if (elem_num == 0)
		return -1;

	clear_bg_mesh();
	pe_buffer.reset();

	grid_hx = hx;
	grid_hy = hy;
	pe_buffer.set_page_size(elem_num);

	// get bounding box
	grid_x_min = nodes[0].x;
	grid_x_max = grid_x_min;
	grid_y_min = nodes[0].y;
	grid_y_max = grid_y_min;
	for (size_t n_id = 1; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		if (grid_x_min > n.x)
			grid_x_min = n.x;
		if (grid_x_max < n.x)
			grid_x_max = n.x;
		if (grid_y_min > n.y)
			grid_y_min = n.y;
		if (grid_y_max < n.y)
			grid_y_max = n.y;
	}

	// init background grid
	double x_padding, y_padding;
	x_padding = (grid_x_max - grid_x_min) * 0.001;
	y_padding = (grid_y_max - grid_y_min) * 0.001;
	grid_x_min -= x_padding;
	grid_x_max += x_padding;
	grid_y_min -= y_padding;
	grid_y_max += y_padding;
	grid_x_num = size_t(ceil((grid_x_max - grid_x_min) / grid_hx));
	grid_y_num = size_t(ceil((grid_y_max - grid_y_min) / grid_hy));
	grid_num = grid_x_num * grid_y_num;
	grid_x_max = grid_x_min + double(grid_x_num) * grid_hx;
	grid_y_max = grid_y_min + double(grid_y_num) * grid_hy;
	bg_grids = new Grid[grid_num];
	size_t grid_id = 0;
	for (size_t row_id = 0; row_id < grid_x_num; ++row_id)
		for (size_t col_id = 0; col_id < grid_y_num; ++col_id)
		{
			Grid &g = bg_grids[grid_id];
			g.x_id = col_id;
			g.y_id = row_id;
			g.pelems = nullptr;
			++grid_id;
		}

	// add element to grids
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
		add_elem_to_bg_grid(elems[e_id]);
}

inline void swap(double &a, double &b)
{
	double c = a;
	a = b;
	b = c;
}

bool test_AABB_triangle_intersection(double xl, double xu, double yl, double yu,
	double x0, double y0, double x1, double y1, double x2, double y2)
{
	double hx, hy, xc, yc;
	hx = xu - xl;
	hy = yu - yl;
	xc = (xl + xu) * 0.5;
	yc = (yl + yu) * 0.5;
	// translate triangle
	// take grid centre as origin
	x0 -= xc;
	y0 -= yc;
	x1 -= xc;
	y1 -= yc;
	x2 -= xc;
	y2 -= yc;

	double r, v_min, v_max;
	// a31
	r = (hx * abs(y1 - y0) + hy * abs(x1 - x0)) * 0.5;
	v_min = x1 * y0 - x0 * y1;
	v_max = (y0 - y1) * x2 + (x1 - x0) * y2;
	if (v_min > v_max)
		swap(v_min, v_max);
	if (v_min > r || v_max < -r)
		return false;
	// a32
	r = (hx * abs(y2 - y1) + hy * abs(x2 - x1)) * 0.5;
	v_min = x2 * y1 - x1 * y2;
	v_max = (y1 - y2) * x0 + (x2 - x1) * y0;
	if (v_min > v_max)
		swap(v_min, v_max);
	if (v_min > r || v_max < -r)
		return false;
	// a33
	r = (hx * abs(y0 - y2) + hy * abs(x0 - x2)) * 0.5;
	v_min = x0 * y2 - x2 * y0;
	v_max = (y2 - y0) * x1 + (x0 - x2) * y1;
	if (v_min > v_max)
		swap(v_min, v_max);
	if (v_min > r || v_max < -r)
		return false;

	return true;
}

void Model_T2D_CHM_s::add_elem_to_bg_grid(Element &e)
{
	// get bounding box of element
	double e_x_min, e_x_max, e_y_min, e_y_max;
	// node 1
	Node &n1 = nodes[e.n1];
	e_x_min = n1.x;
	e_x_max = e_x_min;
	e_y_min = n1.y;
	e_y_max = e_y_min;
	// node 2
	Node &n2 = nodes[e.n2];
	if (e_x_min > n2.x)
		e_x_min = n2.x;
	if (e_x_max < n2.x)
		e_x_max = n2.x;
	if (e_y_min > n2.y)
		e_y_min = n2.y;
	if (e_y_max < n2.y)
		e_y_max = n2.y;
	// node 3
	Node &n3 = nodes[e.n3];
	if (e_x_min > n3.x)
		e_x_min = n3.x;
	if (e_x_max < n3.x)
		e_x_max = n3.x;
	if (e_y_min > n3.y)
		e_y_min = n3.y;
	if (e_y_max < n3.y)
		e_y_max = n3.y;

	// test elem intersection with bg grid
	size_t min_x_id, max_x_id, min_y_id, max_y_id, grid_row_id;
	double grid_xl, grid_xu, grid_yl, grid_yu;
	min_x_id = size_t(floor((e_x_min - grid_x_min) / grid_hx));
	max_x_id = size_t( ceil((e_x_max - grid_x_min) / grid_hx)) + 1;
	min_y_id = size_t(floor((e_y_min - grid_y_min) / grid_hy));
	max_y_id = size_t( ceil((e_y_max - grid_y_min) / grid_hy)) + 1;
	for (size_t y_id = min_y_id; y_id < max_y_id; ++y_id)
	{
		grid_row_id = grid_x_num * y_id;
		grid_yl = grid_y_min + double(y_id) * grid_hy;
		grid_yu = grid_yl + grid_hy;
		for (size_t x_id = min_x_id; x_id < max_x_id; ++x_id)
		{
			Grid &g = bg_grids[grid_row_id + x_id];
			grid_xl = grid_x_min + double(x_id) * grid_hx;
			grid_xu = grid_xl + grid_hx;
			if (test_AABB_triangle_intersection(grid_xl, grid_xu, grid_yl, grid_yu,
												n1.x, n1.y, n2.x, n2.y, n3.x, n3.y))
				add_elem_to_grid(g, e);
		}
	}
}
