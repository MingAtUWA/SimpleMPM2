#include "SimulationCore_pcp.h"

#include "ParticleContactDetectionGrid.h"

void ParticleContactDetectionGrid::init_grid(
	double xl, double xu,
	double yl, double yu,
	double hx, double hy
	)
{
	if (xl >= xu || yl >= yu)
		return;

	clear_grid();
	x_num = size_t(ceil((xu - xl) / hx));
	y_num = size_t(ceil((yu - yl) / hy));
	x0 = xl;
	xn = xu;
	y0 = yl;
	yn = yu;
	dx = (xu - xl) / double(x_num);
	dy = (yu - yl) / double(y_num);

	grids = new Grid[x_num * y_num];
	Grid *pg = grids;
	for (size_t y_id = 0; y_id < y_num; ++y_id)
		for (size_t x_id = 0; x_id < x_num; ++x_id)
		{
			pg->x_id = x_id;
			pg->y_id = y_id;
			pg->bbox_list = nullptr;
			++pg;
		}

	bbox_mem.reset();
}

inline void ParticleContactDetectionGrid::clear_grid(void)
{
	if (grids)
	{
		delete[] grids;
		grids = nullptr;
	}
	x_num = 0;
	y_num = 0;
}

inline void ParticleContactDetectionGrid::reset_grid(void)
{
	bbox_mem.reset();
	Grid *pg = grids;
	for (size_t y_id = 0; y_id < y_num; ++y_id)
		for (size_t x_id = 0; x_id < x_num; ++x_id)
		{
			pg->bbox_list = nullptr;
			++pg;
		}
}

// return the number of grids that this particle overlaps
size_t ParticleContactDetectionGrid::add_particle(double x, double y, double vol)
{
	BoundingBox bbox(x, y, vol);

	// check whether the particle overlap the mesh
	if (!detect_contact(bbox))
		return 0;

	// trim particle if it lies at the edge
	if (bbox.xl < x0)
		bbox.xl = x0;
	if (bbox.xu > xn)
		bbox.xu = xn;
	if (bbox.yl < y0)
		bbox.yl = y0;
	if (bbox.yu > yn)
		bbox.yu = yn;

	size_t xl_id = size_t(floor((bbox.xl - x0) / dx));
	size_t xu_id = size_t(ceil((bbox.xu - x0) / dx));
	size_t x_num = xu_id - xl_id;
	size_t yl_id = size_t(floor((bbox.yl - y0) / dy));
	size_t yu_id = size_t(ceil((bbox.yu - y0) / dy));
	size_t y_num = yu_id - yl_id;

	Grid *pg;
	if (x_num == 1)
	{
		if (y_num == 1)
		{
			pg = pgrid_by_id(xl_id, yl_id);
			add_bbox_to_grid(bbox, *pg);
			return 1;
		}
		else if (y_num == 2)
		{
			pg = pgrid_by_id(xl_id, yl_id);
			add_bbox_to_grid(bbox, *pg);
			pg += x_num;
			add_bbox_to_grid(bbox, *pg);
			return 2;
		}
	}
	else if (x_num == 2)
	{
		if (y_num == 1)
		{
			pg = pgrid_by_id(xl_id, yl_id);
			add_bbox_to_grid(bbox, *pg);
			++pg;
			add_bbox_to_grid(bbox, *pg);
			return 2;
		}
		else if (y_num == 2)
		{
			pg = pgrid_by_id(xl_id, yl_id);
			add_bbox_to_grid(bbox, *pg);
			++pg;
			add_bbox_to_grid(bbox, *pg);
			pg += x_num;
			add_bbox_to_grid(bbox, *pg);
			--pg;
			add_bbox_to_grid(bbox, *pg);
			return 4;
		}
	}

	size_t elem_num = x_num * y_num;
	for (size_t j = 0; j < y_num; ++j)
	{
		pg = pgrid_by_id(xl_id, yl_id + j);
		for (size_t i = 0; i < x_num; ++i)
		{
			add_bbox_to_grid(bbox, *pg);
			++pg;
		}
	}
	return elem_num;
}

// return value:
// true: in contact
// false: not in contact
bool ParticleContactDetectionGrid::detect_contact(double x, double y, double vol)
{
	BoundingBox bbox(x, y, vol);

	// check whether the particle overlap the mesh
	if (!detect_contact(bbox))
		return 0;

	// trim particle if it lies at the edge
	if (bbox.xl < x0)
		bbox.xl = x0;
	if (bbox.xu > xn)
		bbox.xu = xn;
	if (bbox.yl < y0)
		bbox.yl = y0;
	if (bbox.yu > yn)
		bbox.yu = yn;

	size_t xl_id = size_t(floor((bbox.xl - x0) / dx));
	size_t xu_id = size_t(ceil((bbox.xu - x0) / dx));
	size_t x_num = xu_id - xl_id;
	size_t yl_id = size_t(floor((bbox.yl - y0) / dy));
	size_t yu_id = size_t(ceil((bbox.yu - y0) / dy));
	size_t y_num = yu_id - yl_id;

	Grid *pg;
	if (x_num == 1)
	{
		if (y_num == 1)
		{
			pg = pgrid_by_id(xl_id, yl_id);
			return pg->detect_contact(bbox);
		}
		else if (y_num == 2)
		{
			pg = pgrid_by_id(xl_id, yl_id);
			if (pg->detect_contact(bbox)) return true;
			pg += x_num;
			return pg->detect_contact(bbox);
		}
	}
	else if (x_num == 2)
	{
		if (y_num == 1)
		{
			pg = pgrid_by_id(xl_id, yl_id);
			if (pg->detect_contact(bbox)) return true;
			++pg;
			return pg->detect_contact(bbox);
		}
		else if (y_num == 2)
		{
			pg = pgrid_by_id(xl_id, yl_id);
			if (pg->detect_contact(bbox)) return true;
			++pg;
			if (pg->detect_contact(bbox)) return true;
			pg += x_num;
			if (pg->detect_contact(bbox)) return true;
			--pg;
			return pg->detect_contact(bbox);
		}
	}

	size_t elem_num = x_num * y_num;
	for (size_t j = 0; j < y_num; ++j)
	{
		pg = pgrid_by_id(xl_id, yl_id + j);
		for (size_t i = 0; i < x_num; ++i)
		{
			if (pg->detect_contact(bbox))
				return true;
			++pg;
		}
	}
	return false;
}