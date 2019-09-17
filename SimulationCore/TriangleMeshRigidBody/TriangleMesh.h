#ifndef __TRIANGLE_MESH_H__
#define __TRIANGLE_MESH_H__

#include "ItemArray.hpp"
#include "ItemBuffer.hpp"

#include "Geometry.h"

//#define __DEBUG_TRIANGLE_MESH__

template <typename Item>
struct ItemPointer
{
	Item *item;
	ItemPointer *next;
	inline Item &operator *(void) { return *item; }
};

class RigidBody;
class DrawTriangleMesh;

struct TriangleMesh
{
	friend RigidBody;
	friend DrawTriangleMesh;
public:
	struct Node { size_t id; double x, y; };
	struct Element
	{
		size_t id;
		size_t n1, n2, n3;
		double area; // actually is 2 * area, has signed
	};
	struct Edge { size_t n1, n2; };

	class BgGrid
	{
		friend TriangleMesh;

	public:
		enum class PosType
		{
			Outside = 0,
			AtBoundary = 1,
			Inside = 2
		};

		typedef ItemPointer<Element> PElement;
		typedef ItemPointer<Edge> PEdge;

		struct Grid
		{
			size_t x_id, y_id;
			PosType pos_type;

			PElement *elems; // elements covering this grid
			PEdge *edges; // edges covering this grid

			union // coordinates of the centre
			{
				struct { double x_mid, y_mid; };
				Point centre;
			};
			Grid *prev_by_quadtree; // used by BoundaryGirdQuadTree
		};

		double h, xn, yn; // left and lower bound coord = 0.0
		long long x_num, y_num, num;
		long long x_max_id, y_max_id;
		Grid *grids;

		MemoryUtilities::ItemBuffer<PElement> pelem_mem;
		MemoryUtilities::ItemBuffer<PEdge> pedge_mem;

	public:
		typedef MemoryUtilities::ItemArray<long long, 2, 10> IndexArray;

		BgGrid() : grids(nullptr), x_num(0), y_num(0), num(0) {}
		~BgGrid() { clear_grid(); }

		void alloc_grid(long long elem_x_num, long long elem_y_num, double grid_size);
		void clear_grid(void);

		inline Grid &get_grid(size_t y_id, size_t x_id) const { return grids[y_id * x_num + x_id]; }

		// assume the whole line lies within the background mesh
		void get_intersect_points(double x1, double y1, double x2, double y2,
			long long &y_id0, long long &y_idn, IndexArray &x_ids_mem);

		void add_element(Element &elem, Grid &grid);
		void add_edge(Edge &edge, Grid &grid);
		void reset_element_and_edge(void);
	};

protected:
	size_t node_num;
	Node *nodes;
	size_t elem_num;
	Element *elems;

	size_t boundary_edge_num;
	Edge *boundary_edges;

	Rect bounding_box;
	// mass centre at origin
	Rect moved_bounding_box;

	BgGrid bg_grid;

	// geometry properties
	// init in load_mesh()
	double x_mc, y_mc;
	double r_max; // maximum distance from boundary to mass centre
	double area, moi_area; // mass, moment of inertia

public:
	TriangleMesh() :
		nodes(nullptr), node_num(0),
		elems(nullptr), elem_num(0),
		boundary_edges(nullptr), boundary_edge_num(0) {}
	~TriangleMesh()
	{
		clear_nodes();
		clear_elements();
		clear_edges_at_boundary();
	}
	void clear_nodes(void)
	{
		if (nodes) delete[] nodes;
		nodes = nullptr;
		node_num = 0;
	}
	void alloc_nodes(size_t num)
	{
		clear_nodes();
		if (num == 0) return;
		nodes = new Node[num];
		node_num = num;
	}
	void clear_elements(void)
	{
		if (elems) delete[] elems;
		elems = nullptr;
		elem_num = 0;
	}
	void alloc_elements(size_t num)
	{
		clear_elements();
		if (num == 0) return;
		elems = new Element[num];
		elem_num = num;
	}
	void clear_edges_at_boundary(void)
	{
		if (boundary_edges) delete[] boundary_edges;
		boundary_edges = nullptr;
		boundary_edge_num = 0;
	}

	inline size_t get_node_num(void) const noexcept { return node_num; }
	inline const Node *get_nodes(void) const noexcept { return nodes; }
	inline size_t get_elem_num(void) const noexcept { return elem_num; }
	inline const Element *get_elems(void) const noexcept { return elems; }
	inline BgGrid &get_bg_grid(void) { return bg_grid; }
	inline double get_x_mc(void) const noexcept { return x_mc; }
	inline double get_y_mc(void) const noexcept { return y_mc; }
	inline double get_(void) const noexcept { return area; }
	inline double get_moi_area(void) const noexcept { return moi_area; }

	bool is_in_triangle(Element &elem, Point &p);
	inline Rect get_display_range(void) const { return { 0.0, bg_grid.xn, 0.0, bg_grid.yn }; }

	int load_mesh(const char *file_name);
	// after load_mesh()
	int get_bounding_box(void);
	// after load_mesh()
	int find_edges_at_boundary(void);
	// after
	//   get_bounding_box()
	//   find_edges_at_boundary()
	int init_bg_grid(double _grid_size = -1.0);

public:
	// get distance to object
	// fail if return false
	// only fail when edge deteriorate into point and p locates exactly at that point
	bool distance_to_edge(Edge &edge, Point &p, double &dist, double &n1, double &n2);

protected:
	struct RNodePair
	{
		double r; long long x_id, y_id;
		inline void swap(RNodePair &another)
		{
			double r_tmp = r;
			r = another.r;
			another.r = r_tmp;
			long long id_tmp;
			id_tmp = x_id;
			x_id = another.x_id;
			another.x_id = id_tmp;
			id_tmp = y_id;
			y_id = another.y_id;
			another.y_id = id_tmp;
		}
	};
	size_t tree_height; // init in init_bg_grid()
	Rect bg_grid_rect;

	// helper functions
	inline void init_tree_height(size_t elem_x_num, size_t elem_y_num)
	{
		tree_height = 0;
		size_t tree_width = elem_x_num > elem_y_num ? elem_x_num : elem_y_num;
		if (tree_width == 0)
			return;
		--tree_width;
		while (tree_width)
		{
			tree_width = tree_width >> 1;
			++tree_height;
		}
	}
	// searth the background grid as an quad-tree
	// in the future, build an quad-tree from this background grid may be an good idea
	// return false if no boundary line is found
	bool search_node(Point &p, long long raw_x_id, long long raw_y_id, size_t height,
					 double &dist, double &nx, double &ny, double dist_max);

public: // get distance to the mesh boundary
			// return value: 
			//	1 - outside mesh
			//	0 - inside mesh
			// -1 - no boundary lines (within dist_max) is found
	int distance_to_boundary(Point &p, double &dist, double &nx, double &ny, double dist_max = std::numeric_limits<double>::max());

#ifdef __DEBUG_TRIANGLE_MESH__
public: // for debug purpose
	Edge *closest_edge; // for debug purpose
	void dump_mesh_info(const char *dump_file_name);
	inline void print_node(long long x_id, long long y_id,
						   size_t height, bool change_line = true)
	{
		size_t stride = 1ULL << height;
		std::cout << std::string(tree_height - height, ' ')
				  << "node (" << x_id << ", " << x_id + stride << ", "
				  << y_id << ", " << y_id + stride << ") depth: "
				  << tree_height - height;
		if (change_line) std::cout << "\n";
	}
#endif
};

inline double distance(TriangleMesh::Node &n1, TriangleMesh::Node &n2) noexcept
{
	double x_diff = n1.x - n2.x;
	double y_diff = n1.y - n2.y;
	return sqrt(x_diff * x_diff + y_diff * y_diff);
}

#endif