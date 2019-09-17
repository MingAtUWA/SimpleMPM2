#include "SimulationCore_pcp.h"

#include "TriangleMesh.h"
//#include "BgGrid.h"

void TriangleMesh::BgGrid::alloc_grid(long long elem_x_num,
						long long elem_y_num, 
						double grid_size)
{
	x_num = elem_x_num;
	y_num = elem_y_num;
	num = x_num * y_num;
	x_max_id = x_num - 1;
	y_max_id = y_num - 1;
	h = grid_size;
	xn = double(x_num) * h;
	yn = double(y_num) * h;
	if (grids) delete[] grids;
	grids = new Grid[num];
	size_t k = 0;
	for (long long j = 0; j < y_num; ++j)
		for (long long i = 0; i < x_num; ++i)
		{
			Grid &grid = grids[k];
			grid.x_id = i;
			grid.y_id = j;
			grid.pos_type = PosType::Outside;
			grid.elems = nullptr;
			grid.edges = nullptr;
			grid.x_mid = (double(i) + 0.5) * h;
			grid.y_mid = (double(j) + 0.5) * h;
			++k;
		}
	pelem_mem.reset();
	pedge_mem.reset();
}

void TriangleMesh::BgGrid::clear_grid(void)
{
	if (grids)
	{
		delete[] grids;
		grids = nullptr;
	}
	x_num = 0;
	y_num = 0;
	num = 0;
}

void TriangleMesh::BgGrid::get_intersect_points(
	double x1, double y1, double x2, double y2,
	long long &y_id0, long long &y_idn, IndexArray &x_ids_mem)
{
	long long *x_ids;
	x_ids_mem.reset();

	if (y1 == y2) // horizontal line
	{
		y_id0 = long long(y1 / h);
		y_idn = y_id0 + 1;
		x_ids = x_ids_mem.alloc(2);
		x_ids[0] = long long(x1 / h);
		x_ids[1] = long long(x2 / h);
		return;
	}

	// ensure that y1 < y2
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

	y_id0 = long long(y1 / h);
	y_idn = long long(y2 / h);
	if (double(y_idn) * h < y2)
		++y_idn;
	
	size_t x_id_num = y_idn - y_id0;
	double xi, yi, k = (x2 - x1) / (y2 - y1);
	x_ids_mem.reserve(x_id_num + 1);
	long long x_id;
	x_id = long long(x1 / h);
	x_ids_mem.add(x_id);
	for (size_t i = 1; i < x_id_num; ++i)
	{
		yi = double(y_id0 + i) * h;
		xi = x1 + k * (yi - y1);
		x_id = long long(xi / h);
		x_ids_mem.add(x_id);
	}
	x_id = long long(x2 / h);
	x_ids_mem.add(x_id);
	return;
}

void TriangleMesh::BgGrid::add_element(TriangleMesh::Element &elem, Grid &grid)
{
	PElement *pelem = pelem_mem.alloc();
	pelem->item = &elem;
	pelem->next = grid.elems;
	grid.elems = pelem;
}

void TriangleMesh::BgGrid::add_edge(TriangleMesh::Edge &edge, Grid &grid)
{
	PEdge *pedge = pedge_mem.alloc();
	pedge->item = &edge;
	pedge->next = grid.edges;
	grid.edges = pedge;
}

void TriangleMesh::BgGrid::reset_element_and_edge(void)
{
	pelem_mem.reset();
	pedge_mem.reset();
	for (long long grid_id = 0; grid_id < num; ++grid_id)
	{
		grids[grid_id].elems = nullptr;
		grids[grid_id].edges = nullptr;
	}
}
