#include "SimulationCore_pcp.h"

#include "Geometry.h"
#include "Model_S2D_ME_s_RigidBody.h"

Model_S2D_ME_s_RigidBody::Model_S2D_ME_s_RigidBody() :
	Model("Model_S2D_ME_MPM_s"),
	nodes(nullptr), node_x_num(0), node_y_num(0),
	elems(nullptr), elem_x_num(0), elem_y_num(0),
	pcls(nullptr), pcl_num(0),
	bfx_num(0), bfy_num(0), bfxs(nullptr), bfys(nullptr),
	tx_num(0), ty_num(0), txs(nullptr), tys(nullptr),
	ax_num(0), ay_num(0), axs(nullptr), ays(nullptr),
	vx_num(0), vy_num(0), vxs(nullptr), vys(nullptr) {}

Model_S2D_ME_s_RigidBody::~Model_S2D_ME_s_RigidBody()
{
	clear_mesh();
	clear_pcl();
	if (bfxs) delete[] bfxs;
	bfxs = nullptr;
	bfx_num = 0;
	if (bfys) delete[] bfys;
	bfys = nullptr;
	bfy_num = 0;
	if (txs) delete[] txs;
	txs = nullptr;
	tx_num = 0;
	if (tys) delete[] tys;
	tys = nullptr;
	ty_num = 0;
	if (axs) delete[] axs;
	axs = nullptr;
	ax_num = 0;
	if (ays) delete[] ays;
	ays = nullptr;
	ay_num = 0;
	if (vxs) delete[] vxs;
	vxs = nullptr;
	vx_num = 0;
	if (vys) delete[] vys;
	vys = nullptr;
	vy_num = 0;
}

void Model_S2D_ME_s_RigidBody::init_mesh(
	double grid_size, size_t _elem_x_num, size_t _elem_y_num,
	double x_start, double y_start)
{
	clear_mesh();
	h = grid_size;
	x0 = x_start;
	xn = x0 + h * double(_elem_x_num);
	y0 = y_start;
	yn = y0 + h * double(_elem_y_num);
	size_t i, j, k;
	// elements
	elem_x_num = _elem_x_num;
	elem_y_num = _elem_y_num;
	elem_num = elem_x_num * elem_y_num;
	elems = new Element[elem_num];
	k = 0;
	for (i = 0; i < elem_y_num; ++i)
		for (j = 0; j < elem_x_num; ++j)
		{
			elems[k].index_x = j;
			elems[k].index_y = i;
			++k;
		}
	// nodes
	node_x_num = elem_x_num + 1;
	node_y_num = elem_y_num + 1;
	node_num = node_x_num * node_y_num;
	nodes = new Node[node_num];
	k = 0;
	for (i = 0; i < node_y_num; ++i)
		for (j = 0; j < node_x_num; ++j)
		{
			nodes[k].index_x = j;
			nodes[k].index_y = i;
			++k;
		}
}

void Model_S2D_ME_s_RigidBody::clear_mesh(void)
{
	if (nodes) delete[] nodes;
	nodes = nullptr;
	node_x_num = 0;
	node_y_num = 0;
	node_num = 0;
	if (elems) delete[] elems;
	elems = nullptr;
	elem_x_num = 0;
	elem_y_num = 0;
	elem_num = 0;
}

void Model_S2D_ME_s_RigidBody::init_pcl(
	size_t num, double m, double density, double E, double niu)
{
	clear_pcl();
	pcl_num = num;
	pcls = new Particle[num];
	for (size_t i = 0; i < num; ++i)
	{
		Particle &pcl = pcls[i];
		pcl.index = i;
		pcl.x = 0.0;
		pcl.y = 0.0;
		pcl.vx = 0.0;
		pcl.vy = 0.0;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
		pcl.m = m;
		pcl.density = density;
		pcl.E = E;
		pcl.niu = niu;
		pcl.s11 = 0.0;
		pcl.s12 = 0.0;
		pcl.s22 = 0.0;
		pcl.e11 = 0.0;
		pcl.e12 = 0.0;
		pcl.e22 = 0.0;
	}
}

void Model_S2D_ME_s_RigidBody::clear_pcl(void)
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
	}
	pcl_num = 0;
}

bool Model_S2D_ME_s_RigidBody::get_intersect_points(
		double x1, double y1, double x2, double y2,
		LongLongRange &y_id_range, PreAllocIdArray &x_ids_mem)
{
	long long *x_ids;
	size_t y_id0, y_idn;
	if (y1 == y2)
	{
		if (y1 < y0 || y1 >= yn)
			return false;
		y_id0 = long long((y1 - y0) / h);
		y_idn = y_id0 + 1;
		if (y_idn > elem_y_num)
			y_idn = elem_y_num;
		y_id_range.lower = y_id0;
		y_id_range.upper = y_idn;
		x_ids_mem.reset();
		x_ids = x_ids_mem.alloc(2);
		x_ids[0] = long long(floor((x1 - x0) / h));
		x_ids[1] = long long(floor((x2 - x0) / h));
		return true;
	}

	// y1 != y2
	if (y1 > y2)
	{
		double d_tmp;
		d_tmp = y1;
		y1 = y2;
		y2 = d_tmp;
		d_tmp = x1;
		x1 = x2;
		x2 = d_tmp;
	}
	if (y2 <= y0 || y1 >= yn)
		return false;
	double xi, yi, k = (x2 - x1) / (y2 - y1);
	size_t x_id0, x_idn;
	if (y1 < y0)
	{
		y_id0 = 0;
		xi = x1 + k * (y0 - y1);
		x_id0 = long long(floor((xi - x0) / h));
	}
	else
	{
		y_id0 = size_t((y1 - y0) / h);
		x_id0 = long long(floor((x1 - x0) / h));
	}
	if (y2 >= yn)
	{
		y_idn = elem_y_num;
		xi = x1 + k * (yn - y1);
		x_idn = long long(floor((xi - x0) / h));
	}
	else
	{
		y_idn = size_t((y2 - y0) / h);
		if (double(y_idn) * h < y2 - y0)
			++y_idn;
		x_idn = long long(floor((x2 - y0) / h));
	}
	y_id_range.lower = y_id0;
	y_id_range.upper = y_idn;
	size_t x_id_num = y_idn - y_id0 + 1;
	x_ids_mem.reset();
	x_ids = x_ids_mem.alloc(x_id_num);
	// the first x_id
	x_ids[0] = x_id0;
	// the last x_id
	x_ids[x_id_num-1] = x_idn;
	// the other x_id
	for (size_t i = 1; i < x_id_num-1; ++i)
	{
		yi = y0 + double(y_id0 + i) * h;
		xi = x1 + k * (yi - y1);
		x_ids[i] = long long(floor((xi - x0) / h));
	}
	return true;
}

void Model_S2D_ME_s_RigidBody::update_x_id_range(
	LongLongRange &y_id_range, PreAllocIdArray & x_ids_mem)
{
	size_t start_id = y_id_range.lower - y_id_min;
	size_t id_num = y_id_range.upper - y_id_range.lower;
	long long x_id_min, x_id_max;
	long long *x_ids = x_ids_mem.get_mem();
	LongLongRange *x_range = x_range_mem.get_mem();
	for (size_t i = start_id; i < id_num; ++i)
	{
		if (x_ids[i] < x_ids[i + 1])
		{
			x_id_min = x_ids[i];
			x_id_max = x_ids[i + 1] + 1;
		}
		else
		{
			x_id_min = x_ids[i + 1];
			x_id_max = x_ids[i] + 1;
		}

		if (x_range[i].lower > x_id_min)
			x_range[i].lower = x_id_min;
		if (x_range[i].upper < x_id_max)
			x_range[i].upper = x_id_max;
	}
}

// clip and resterize, need better algorithm later
void Model_S2D_ME_s_RigidBody::rasterize_rect_on_grid(
		double x1, double y1, double x2, double y2,
		double x3, double y3, double x4, double y4)
{
	bool res1, res2, res3, res4;
	// line1
	res1 = get_intersect_points(x1, y1, x2, y2, y_id_range1, x_ids_mem1);
	// line2
	res2 = get_intersect_points(x2, y2, x3, y3, y_id_range2, x_ids_mem2);
	// line3
	res3 = get_intersect_points(x3, y3, x4, y4, y_id_range3, x_ids_mem3);
	// line4
	res4 = get_intersect_points(x4, y4, x1, y1, y_id_range4, x_ids_mem4);

	y_id_min = elem_y_num;
	long long y_id_max = 0;
	if (res1)
	{
		if (y_id_min > y_id_range1.lower)
			y_id_min = y_id_range1.lower;
		if (y_id_max < y_id_range1.upper)
			y_id_max = y_id_range1.upper;
	}
	if (res2)
	{
		if (y_id_min > y_id_range2.lower)
			y_id_min = y_id_range2.lower;
		if (y_id_max < y_id_range2.upper)
			y_id_max = y_id_range2.upper;
	}
	if (res3)
	{
		if (y_id_min > y_id_range3.lower)
			y_id_min = y_id_range3.lower;
		if (y_id_max < y_id_range3.upper)
			y_id_max = y_id_range3.upper;
	}
	if (res4)
	{
		if (y_id_min > y_id_range4.lower)
			y_id_min = y_id_range4.lower;
		if (y_id_max < y_id_range4.upper)
			y_id_max = y_id_range4.upper;
	}
	if (y_id_min > y_id_max)
		return;

	long long x_range_num = y_id_max - y_id_min;
	x_range_mem.reset();
	LongLongRange *x_range = x_range_mem.alloc(x_range_num);
	// init x_range
	for (long long range_id = 0; range_id < x_range_num; ++range_id)
	{
		x_range[range_id].lower = long long(elem_x_num);
		x_range[range_id].upper = 0;
	}
	// set x_range
	if (res1) update_x_id_range(y_id_range1, x_ids_mem1);
	if (res2) update_x_id_range(y_id_range2, x_ids_mem2);
	if (res3) update_x_id_range(y_id_range3, x_ids_mem3);
	if (res4) update_x_id_range(y_id_range4, x_ids_mem4);
	// clip x_range
	for (long long range_id = 0; range_id < x_range_num; ++range_id)
	{
		if (x_range[range_id].lower < 0)
			x_range[range_id].lower = 0;
		if (x_range[range_id].upper > long long(elem_x_num))
			x_range[range_id].upper = long long(elem_x_num);
	}

	// rasterize the grid
	for (long long range_id = 0; range_id < x_range_num; ++range_id)
	{
		long long y_id = y_id_min + range_id;
		for (long long x_id = x_range[range_id].lower;
			 x_id < x_range[range_id].upper; ++x_id)
			elems[y_id * elem_x_num + x_id].has_rigid_object = true;
	}
}
