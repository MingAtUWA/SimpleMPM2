#include "SimulationCore_pcp.h"

#include <fstream>

#include "TriangleMesh.h"

#include "HashTable.h"
#include "HashTableEdge.h"

template<typename _Type>
inline void swap(_Type &a, _Type &b) noexcept { _Type tmp = a; a = b, b = tmp; }

bool TriangleMesh::is_in_triangle(Element &elem, Point &p)
{
	Node &n1 = nodes[elem.n1];
	Node &n2 = nodes[elem.n2];
	Node &n3 = nodes[elem.n3];
	double a = (n2.x - n1.x) * (p.y - n1.y) - (p.x - n1.x) * (n2.y - n1.y);
	double b = (p.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (p.y - n1.y);
	double c = elem.area - a - b;
	return 0.0 <= a && a <= elem.area
		&& 0.0 <= b && b <= elem.area 
		&& 0.0 <= c && c <= elem.area;
}

void TriangleMesh::compress_node_and_elem_indices(void)
{
	// "compress" nodes and elements index
	HashTable node_id_map(node_num);
	for (size_t i = 0; i < node_num; ++i)
	{
		node_id_map.add_pair(nodes[i].id, i);
		nodes[i].id = i;
	}
	for (size_t i = 0; i < elem_num; ++i)
	{
		size_t n_id;
		Element &elem = elems[i];
		node_id_map.get_pair(elem.n1, n_id);
		elem.n1 = n_id;
		node_id_map.get_pair(elem.n2, n_id);
		elem.n2 = n_id;
		node_id_map.get_pair(elem.n3, n_id);
		elem.n3 = n_id;
		elem.id = i;
	}
}

void TriangleMesh::cal_area_and_reorder_node(void)
{
	// Calculate area of triangle.
	area = 0.0;
	x_mc = 0.0;
	y_mc = 0.0;
	moi_area = 0.0;
	double elem_area;
	for (size_t i = 0; i < elem_num; ++i)
	{
		Element &elem = elems[i];
		Node &n1 = nodes[elem.n1];
		Node &n2 = nodes[elem.n2];
		Node &n3 = nodes[elem.n3];
		double x1, y1, x2, y2, x3, y3;
		x1 = n1.x;
		y1 = n1.y;
		x2 = n2.x;
		y2 = n2.y;
		x3 = n3.x;
		y3 = n3.y;
		elem.area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
		// Ensure that nodes index of element is in counter-clockwise sequence.
		if (elem.area < 0.0)
		{
			elem.area = -elem.area;
			size_t n_tmp = elem.n2;
			elem.n2 = elem.n3;
			elem.n3 = n_tmp;
		}

		elem_area = elem.area / 2.0;
		area += elem_area;
		x_mc += elem_area * (x1 + x2 + x3) / 3.0;
		y_mc += elem_area * (y1 + y2 + y3) / 3.0;
		moi_area += (x1 * x1 + x2 * x2 + x3 * x3
			+ x1 * x2 + x2 * x3 + x3 * x1
			+ y1 * y1 + y2 * y2 + y3 * y3
			+ y1 * y2 + y2 * y3 + y3 * y1) * elem_area / 6.0;
	}
	x_mc /= area;
	y_mc /= area;
	moi_area -= area * (x_mc * x_mc + y_mc * y_mc);
}

int TriangleMesh::load_mesh(const char *file_name)
{
	union
	{
		char num_buffer[8];
		unsigned long long num;
	};

	union
	{
		char node_buffer[32];
		struct
		{
			unsigned long long id;
			double x, y, z;
		} node_info;
	};

	union
	{
		char elem_buffer[32];
		struct { unsigned long long id, n1, n2, n3; } elem_info;
	};

	std::ifstream mesh_file(file_name, std::ios::binary);
	if (!mesh_file.is_open())
		return -1;

	// load node
	mesh_file.read(num_buffer, sizeof(num_buffer));
	alloc_nodes(num);
	for (size_t i = 0; i < num; ++i)
	{
		mesh_file.read(node_buffer, sizeof(node_buffer));
		nodes[i].id = node_info.id;
		nodes[i].x = node_info.x;
		nodes[i].y = node_info.y;
	}

	// load element
	mesh_file.read(num_buffer, sizeof(num_buffer));
	alloc_elements(num);
	for (size_t i = 0; i < num; ++i)
	{
		mesh_file.read(elem_buffer, sizeof(elem_buffer));
		elems[i].id = elem_info.id;
		elems[i].n1 = elem_info.n1;
		elems[i].n2 = elem_info.n2;
		elems[i].n3 = elem_info.n3;
	}

	if (node_num == 0 || elem_num == 0)
		return -2;

	compress_node_and_elem_indices();
	cal_area_and_reorder_node();

	return 0;
}

int TriangleMesh::load_mesh(double *node_coords, size_t node_num,
							size_t *elem_indices, size_t elem_num)
{
	if (node_num == 0 || elem_num == 0)
		return -2;

	size_t max_node_id;
	max_node_id = elem_indices[0];
	for (size_t i = 1; i < 3*elem_num; ++i)
	{
		if (max_node_id < elem_indices[i])
			max_node_id = elem_indices[i];
	}
	if (max_node_id >= node_num)
		return -2;

	alloc_nodes(node_num);
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		nodes[n_id].id = n_id;
		nodes[n_id].x = node_coords[2 * n_id];
		nodes[n_id].y = node_coords[2 * n_id + 1];
	}
	alloc_elements(elem_num);
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		elems[e_id].id = e_id;
		elems[e_id].n1 = elem_indices[3 * e_id];
		elems[e_id].n2 = elem_indices[3 * e_id + 1];
		elems[e_id].n3 = elem_indices[3 * e_id + 2];
	}

	cal_area_and_reorder_node();

	return 0;

}

int TriangleMesh::get_bounding_box(void)
{
	if (node_num == 0)
		return -1;

	// Init bounding box
	bounding_box.xl = nodes[0].x;
	bounding_box.xu = bounding_box.xl;
	bounding_box.yl = nodes[0].y;
	bounding_box.yu = bounding_box.yl;
	for (size_t i = 1; i < node_num; ++i)
	{
		if (bounding_box.xl > nodes[i].x)
			bounding_box.xl = nodes[i].x;
		if (bounding_box.xu < nodes[i].x)
			bounding_box.xu = nodes[i].x;
		if (bounding_box.yl > nodes[i].y)
			bounding_box.yl = nodes[i].y;
		if (bounding_box.yu < nodes[i].y)
			bounding_box.yu = nodes[i].y;
	}

	moved_bounding_box.xl = bounding_box.xl - x_mc;
	moved_bounding_box.xu = bounding_box.xu - x_mc;
	moved_bounding_box.yl = bounding_box.yl - y_mc;
	moved_bounding_box.yu = bounding_box.yu - y_mc;

	return 0;
}

int TriangleMesh::find_edges_at_boundary(void)
{
	clear_edges_at_boundary();
	HashTableEdge hash_edge(node_num);
	for (size_t edge_id = 0; edge_id < elem_num; ++edge_id)
	{
		Element &elem = elems[edge_id];
		hash_edge.add_edge(elem.n1, elem.n2);
		hash_edge.add_edge(elem.n2, elem.n3);
		hash_edge.add_edge(elem.n3, elem.n1);
	}
	boundary_edge_num = hash_edge.get_edge_at_boundary(boundary_edges);
	if (boundary_edge_num == 0)
		return -1;
	return 0;
}

// After edges at boundaries are obtained
int TriangleMesh::init_bg_grid(double _grid_size)
{
	double grid_size;
	if (_grid_size > 0)
	{
		grid_size = _grid_size;
	}
	else
	{
		// Use length of shortest boundary segment as size of background grid
		double edge_len_min, edge_len_tmp;
		size_t be_id;
		for (be_id = 0; be_id < boundary_edge_num; ++be_id)
		{
			Edge &be = boundary_edges[be_id];
			edge_len_min = distance(nodes[be.n1], nodes[be.n2]);
			if (edge_len_min != 0.0)
				break;
		}
		for (; be_id < boundary_edge_num; ++be_id)
		{
			Edge &be = boundary_edges[be_id];
			edge_len_tmp = distance(nodes[be.n1], nodes[be.n2]);
			if (edge_len_tmp != 0.0 && edge_len_min > edge_len_tmp)
				edge_len_min = edge_len_tmp;
		}
		if (edge_len_min == 0)
			return -1;
		grid_size = edge_len_min;
	}

	// Calculate number of element in x and y direction
	double bb_len_x = bounding_box.xu - bounding_box.xl;
	double bb_len_y = bounding_box.yu - bounding_box.yl;
	size_t elem_x_num = size_t(bb_len_x / grid_size) + 1;
	size_t elem_y_num = size_t(bb_len_y / grid_size) + 1;
	// init the tree height
	init_tree_height(elem_x_num, elem_y_num);
	// init grid
	bg_grid.alloc_grid(elem_x_num, elem_y_num, grid_size);
	bg_grid_rect.xl = 0.0;
	bg_grid_rect.xu = bg_grid.xn;
	bg_grid_rect.yl = 0.0;
	bg_grid_rect.yu = bg_grid.yn;

	// Adjust coordinates of bounding box and nodes
	double offset_x, offset_y;
	offset_x = (double(elem_x_num) * grid_size - bb_len_x) * 0.5 - bounding_box.xl;
	offset_y = (double(elem_y_num) * grid_size - bb_len_y) * 0.5 - bounding_box.yl;
	bounding_box.xl += offset_x;
	bounding_box.xu += offset_x;
	bounding_box.yl += offset_y;
	bounding_box.yu += offset_y;
	x_mc += offset_x;
	y_mc += offset_y;
	//std::cout << x_mc << ", " << y_mc << "\n";
	for (size_t i = 0; i < node_num; ++i)
	{
		Node &n = nodes[i];
		n.x += offset_x;
		n.y += offset_y;
	}

	// Get element at boundary
	double x1, y1, x2, y2, x3, y3;
	long long y_id0, y_idn, *x_ids;
	BgGrid::IndexArray x_ids_mem;
	size_t line_num, x_id_min, x_id_max;
	for (size_t edge_id = 0; edge_id < boundary_edge_num; ++edge_id)
	{
		Edge &edge = boundary_edges[edge_id];
		Node &n1 = nodes[edge.n1];
		x1 = n1.x, y1 = n1.y;
		Node &n2 = nodes[edge.n2];
		x2 = n2.x, y2 = n2.y;
		bg_grid.get_intersect_points(x1, y1, x2, y2, y_id0, y_idn, x_ids_mem);
		x_ids = x_ids_mem.get_mem();
		line_num = y_idn - y_id0;
		for (size_t line_id = 0; line_id < line_num; ++line_id)
		{
			if (x_ids[line_id] < x_ids[line_id + 1])
			{
				x_id_min = x_ids[line_id];
				x_id_max = x_ids[line_id + 1];
			}
			else
			{
				x_id_min = x_ids[line_id + 1];
				x_id_max = x_ids[line_id];
			}
			++x_id_max;
			for (size_t col_id = x_id_min; col_id < x_id_max; ++col_id)
			{
				BgGrid::Grid &grid = bg_grid.get_grid(y_id0 + line_id, col_id);
				grid.pos_type = BgGrid::PosType::AtBoundary;
				bg_grid.add_edge(edge, grid);
			}
		}
	}

	// Get element inside
	long long *x_id_range;
	BgGrid::IndexArray x_id_range_mem;
	for (size_t elem_id = 0; elem_id < elem_num; ++elem_id)
	{
		Element &elem = elems[elem_id];
		// find the highest (y1), middle (y2) and lowest (y3) point
		Node &n1 = nodes[elem.n1];
		x1 = n1.x; y1 = n1.y;
		Node &n2 = nodes[elem.n2];
		x2 = n2.x; y2 = n2.y;
		Node &n3 = nodes[elem.n3];
		x3 = n3.x; y3 = n3.y;
		if (y1 < y2)
		{
			swap(y1, y2);
			swap(x1, x2);
		}
		if (y1 < y3)
		{
			swap(y1, y3);
			swap(x1, x3);
		}
		if (y3 > y2)
		{
			swap(y2, y3);
			swap(x2, x3);
		}
		// line 13 (widest span in y direction)
		size_t y13_num, y23_num, y12_num;
		size_t y_start;
		bg_grid.get_intersect_points(x3, y3, x1, y1, y_id0, y_idn, x_ids_mem);
		//for (size_t i = 0; i < x_ids_mem.get_num(); i++)
		//	std::cout << x_ids_mem[i] << ", ";
		//std::cout << "\n";
		y_start = y_id0;
		y13_num = y_idn - y_id0 + 1;
		x_id_range_mem.reset();
		x_id_range = x_id_range_mem.alloc(y13_num + y13_num);
		for (size_t i = 0; i < y13_num; ++i)
		{
			x_id_range[2 * i] = x_ids[i];
			x_id_range[2 * i + 1] = x_ids[i];
		}
		// line 23
		bg_grid.get_intersect_points(x3, y3, x2, y2,
									 y_id0, y_idn, x_ids_mem);
		x_ids = x_ids_mem.get_mem();
		y23_num = y_idn - y_id0 + 1;
		for (size_t i = 0; i < y23_num; ++i)
		{
			if (x_id_range[2*i] > x_ids[i])
				x_id_range[2*i] = x_ids[i];
			if (x_id_range[2*i+1] < x_ids[i])
				x_id_range[2*i+1] = x_ids[i];
		}
		// line 12
		bg_grid.get_intersect_points(x2, y2, x1, y1, y_id0, y_idn, x_ids_mem);
		x_ids = x_ids_mem.get_mem();
		y12_num = y_idn - y_id0 + 1;
		size_t y12_start = y13_num - y12_num;
		for (size_t i = 0; i < y12_num; ++i)
		{
			if (x_id_range[2*(y12_start+i)] > x_ids[i])
				x_id_range[2*(y12_start+i)] = x_ids[i];
			if (x_id_range[2*(y12_start+i)+1] < x_ids[i])
				x_id_range[2*(y12_start+i)+1] = x_ids[i];
		}

		// get grid inside the triangular mesh
		size_t y_num = y13_num - 1;
		size_t x_min, x_max;
		for (size_t line_id = 0; line_id < y_num; ++line_id)
		{
			x_min = x_id_range[line_id*2] < x_id_range[line_id*2+2] ? x_id_range[line_id*2] : x_id_range[line_id*2+2];
			x_max = x_id_range[line_id*2+1] > x_id_range[line_id*2+3] ? x_id_range[line_id*2+1] : x_id_range[line_id*2+3];
			++x_max;
			for (size_t col_id = x_min; col_id < x_max; ++col_id)
			{
				BgGrid::Grid &grid = bg_grid.get_grid(y_start + line_id, col_id);
				if (grid.pos_type == BgGrid::PosType::AtBoundary)
				{
					bg_grid.add_element(elem, grid);
				}
				else
				{
					grid.pos_type = BgGrid::PosType::Inside;
				}
			}
		}
	}

	return 0;
}

inline bool is_lie_within(double x1, double x2, double x)
{
	return x1 > x2 ? x2 <= x && x <= x1 : x1 <= x && x <= x2;
}

bool TriangleMesh::distance_to_edge(Edge &edge, Point &p, 
									double &dist, double &nx, double &ny)
{
	Node &n1 = nodes[edge.n1];
	Node &n2 = nodes[edge.n2];
	double x1, y1, x2, y2, x, y;
	x1 = n1.x;
	y1 = n1.y;
	x2 = n2.x;
	y2 = n2.y;
	x = p.x;
	y = p.y;

	double x_diff, y_diff, len;
	if (x1 == x2 && y1 == y2) // the line deteriorates into a point
	{
		x_diff = x - x1;
		y_diff = y - y1;
		dist = sqrt(x_diff * x_diff + y_diff * y_diff);
		if (dist == 0)
			return false;
		nx = x_diff / dist;
		ny = y_diff / dist;
		return true;
	}

	x_diff = x2 - x1;
	y_diff = y2 - y1;
	len = sqrt(x_diff * x_diff + y_diff * y_diff);
	nx = -y_diff / len;
	ny = x_diff / len;
	dist = ((y - y2) * (x1 - x2) - (x - x2) * (y1 - y2)) / len;
	double x_interset = x + nx * dist;
	double y_interset = y + ny * dist;
	if (dist > 0)
	{
		nx = -nx;
		ny = -ny;
	}
	else
	{
		dist = -dist;
	}
	// projection point locates within line
	if (is_lie_within(x1, x2, x_interset) &&
		is_lie_within(y1, y2, y_interset))
		return true;
	
	x_diff = x - x1;
	y_diff = y - y1;
	dist = x_diff * x_diff + y_diff * y_diff;
	double x_diff2, y_diff2, dist2;
	x_diff2 = x - x2;
	y_diff2 = y - y2;
	dist2 = x_diff2 * x_diff2 + y_diff2 * y_diff2;
	if (dist < dist2)
	{
		dist = sqrt(dist);
		nx = x_diff / dist;
		ny = y_diff / dist;
	}
	else
	{
		dist = sqrt(dist2);
		nx = x_diff2 / dist;
		ny = y_diff2 / dist;
	}

	return true;
}



int TriangleMesh::distance_to_boundary(Point &p, double &dist, double &nx, double &ny, double dist_max)
{
#ifdef __DEBUG_TRIANGLE_MESH__
	std::cout << "Point (" << p.x << ", " << p.y << ") is ";
#endif

	int res = 1;
	if (p.x < 0.0 || p.x >= bg_grid.xn || p.y < 0.0 || p.y >= bg_grid.yn)
	{
#ifdef __DEBUG_TRIANGLE_MESH__
		std::cout << "outside mesh.\n";
#endif
		if (distance(bg_grid_rect, p) < dist_max)
			goto find_closest_boundary_grid;
		return -1;
	}

	BgGrid::Grid &in_grid = bg_grid.get_grid(size_t(p.y / bg_grid.h), size_t(p.x / bg_grid.h));
	switch (in_grid.pos_type)
	{
	case BgGrid::PosType::Outside:
#ifdef __DEBUG_TRIANGLE_MESH__
		std::cout << "outside mesh.\n";
#endif
		break;
	case BgGrid::PosType::AtBoundary:
#ifdef __DEBUG_TRIANGLE_MESH__
		std::cout << "close to mesh boundary.\n";
#endif
		for (BgGrid::PElement *pelem = in_grid.elems;
			pelem; pelem = pelem->next)
		{
			if (is_in_triangle(*(pelem->item), p))
			{
				res = 0;
				break;
			}
		}
		break;
	case BgGrid::PosType::Inside:
#ifdef __DEBUG_TRIANGLE_MESH__
		std::cout << "inside mesh.\n";
#endif
		res = 0;
		break;
	default:
#ifdef __DEBUG_TRIANGLE_MESH__
		// unknow grid type
		throw std::exception("\nError: unknown background grid position\n");
#endif
		return -2;
	}
	
find_closest_boundary_grid:
	dist = dist_max;
#ifdef __DEBUG_TRIANGLE_MESH__
	std::cout << "Start searching...\n";
#endif
	if (search_node(p, 0, 0, tree_height, dist, nx, ny, dist_max))
	{
		if (res == 0) // inside 
		{
			nx = -nx;
			ny = -ny;
		}
		else // outide
		{
			dist = -dist;
		}

#ifdef __DEBUG_TRIANGLE_MESH__
		Node &ce_n1 = nodes[closest_edge->n1], &ce_n2 = nodes[closest_edge->n2];
		std::cout << "End searching.\nClosest edge is (" << closest_edge->n1 << ", "
				  << closest_edge->n2 << "), from (" << ce_n1.x <<  ", " << ce_n1.y
				  << ") to (" << ce_n2.x << ", " << ce_n2.y << ").\n"
				  "Distance: " << dist << ", normal: (" << nx << ", " << ny << ")\n";;
#endif

		return res;
	}

#ifdef __DEBUG_TRIANGLE_MESH__
	std::cout << "End searching.\nClosest edge not found.\n";
#endif
	return -1;
}


// helper function
// search background grid as a quadtree
bool TriangleMesh::search_node(Point &p, long long raw_x_id0, long long raw_y_id0, size_t height,
							   double &dist, double &nx, double &ny, double dist_max)
{
	if (raw_x_id0 > bg_grid.x_max_id || raw_y_id0 > bg_grid.y_max_id)
	{
#ifdef __DEBUG_TRIANGLE_MESH__
		print_node(raw_x_id0, raw_y_id0, height, false);
		std::cout << " - not in range.\n";
#endif
		return false;
	}

#ifdef __DEBUG_TRIANGLE_MESH__
	print_node(raw_x_id0, raw_y_id0, height);
#endif

	if (height) // non leaf node
	{
		--height; // height is child node
		long long stride, raw_x_id1, raw_y_id1;
		double stride_len, x_low, x_mid, x_up, y_low, y_mid, y_up;
		stride = 1ULL << height;
		raw_x_id1 = raw_x_id0 + stride;
		raw_y_id1 = raw_y_id0 + stride;
		stride_len = stride * bg_grid.h;
		x_low = raw_x_id0 * bg_grid.h;
		x_mid = x_low + stride_len;
		x_up  = x_mid + stride_len;
		y_low = raw_y_id0 * bg_grid.h;
		y_mid = y_low + stride_len;
		y_up  = y_mid + stride_len;
		Rect child_rect;
		// bottom left
		child_rect = { x_low, x_mid, y_low, y_mid };
		RNodePair rn_pair1 = { distance(child_rect, p), raw_x_id0, raw_y_id0 };
		// bottom right
		child_rect = { x_mid,  x_up, y_low, y_mid };
		RNodePair rn_pair2 = { distance(child_rect, p), raw_x_id1, raw_y_id0 };
		// upper left
		child_rect = { x_low, x_mid, y_mid, y_up };
		RNodePair rn_pair3 = { distance(child_rect, p), raw_x_id0, raw_y_id1 };
		// upper right
		child_rect = { x_mid,  x_up, y_mid, y_up };
		RNodePair rn_pair4 = { distance(child_rect, p), raw_x_id1, raw_y_id1 };

		if (rn_pair1.r > rn_pair2.r)
			rn_pair1.swap(rn_pair2);
		if (rn_pair1.r > rn_pair3.r)
			rn_pair1.swap(rn_pair3);
		if (rn_pair1.r > rn_pair4.r)
			rn_pair1.swap(rn_pair4);
		if (rn_pair1.r > dist || rn_pair1.r > dist_max)
			return false;
		bool res = search_node(p, rn_pair1.x_id, rn_pair1.y_id,
							   height, dist, nx, ny, dist_max);

		if (rn_pair2.r > rn_pair3.r)
			rn_pair2.swap(rn_pair3);
		if (rn_pair2.r > rn_pair4.r)
			rn_pair2.swap(rn_pair4);
		if (rn_pair2.r > dist || rn_pair2.r > dist_max)
			return res;
		res = search_node(p, rn_pair2.x_id, rn_pair2.y_id,
						  height, dist, nx, ny, dist_max) || res;

		if (rn_pair3.r > rn_pair4.r)
			rn_pair3.swap(rn_pair4);
		if (rn_pair3.r > dist || rn_pair3.r > dist_max)
			return res;
		res = search_node(p, rn_pair3.x_id, rn_pair3.y_id,
						  height, dist, nx, ny, dist_max) || res;

		if (rn_pair4.r > dist || rn_pair4.r > dist_max)
			return res;
		res = search_node(p, rn_pair4.x_id, rn_pair4.y_id,
						  height, dist, nx, ny, dist_max) || res;

		return res;
	}

	// leaf node
	BgGrid::Grid &cur_grid = bg_grid.get_grid(raw_y_id0, raw_x_id0);
	double dist_tmp, nx_tmp, ny_tmp;
	bool res = false;
	if (cur_grid.pos_type == BgGrid::PosType::AtBoundary)
	{
		for (BgGrid::PEdge *pedge = cur_grid.edges; pedge; pedge = pedge->next)
		{
			if (distance_to_edge(*(pedge->item), p, dist_tmp, nx_tmp, ny_tmp) &&
				dist > dist_tmp)
			{
				dist = dist_tmp;
				nx = nx_tmp;
				ny = ny_tmp;
				res = true;
#ifdef __DEBUG_TRIANGLE_MESH__
				closest_edge = pedge->item;
#endif
			}
		}
	}
	return res;
}

#ifdef __DEBUG_TRIANGLE_MESH__
void TriangleMesh::dump_mesh_info(const char *dump_file_name)
{
	std::fstream dump_file;
	dump_file.open(dump_file_name, std::ios::out | std::ios::binary);

	dump_file << "$Node - id, x, y\n";
	for (size_t i = 0; i < node_num; ++i)
	{
		Node &n = nodes[i];
		dump_file << n.id << ", " << n.x << ", " << n.y << "\n";
	}

	dump_file << "$Element - id, n1, n2, n3\n";
	for (size_t i = 0; i < elem_num; ++i)
	{
		Element &elem = elems[i];
		dump_file << elem.id << ", " << elem.n1 << ", "
			<< elem.n2 << ", " << elem.n3 << "\n";
	}

	dump_file << "$Boundary line - id, n1, n2\n";
	for (size_t i = 0; i < boundary_edge_num; ++i)
	{
		Edge &edge = boundary_edges[i];
		dump_file << i << ", " << edge.n1 << ", " << edge.n2 << "\n";
	}

	dump_file.close();
}
#endif