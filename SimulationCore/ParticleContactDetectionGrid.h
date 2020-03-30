#ifndef __Particle_Contact_Detection_Grid_HPP_
#define __Particle_Contact_Detection_Grid_HPP_

#include <cmath>
#include "ItemBuffer.hpp"

class ParticleContactDetectionGrid
{
protected:
	struct BoundingBox
	{
		double xl, xu, yl, yu;
		BoundingBox(double x, double y, double vol) { init(x, y, vol); }
		inline void init(double x, double y, double vol)
		{
			double pcl_hlen = sqrt(vol) * 0.5;
			xl = x - pcl_hlen;
			xu = x + pcl_hlen;
			yl = y - pcl_hlen;
			yu = y + pcl_hlen;
		}
		inline void copy(BoundingBox &bbox)
		{
			xl = bbox.xl;
			xu = bbox.xu;
			yl = bbox.yl;
			yu = bbox.yu;
		}
		inline bool detect_contact(BoundingBox &bbox)
		{
			return xu <= bbox.xl || xl >= bbox.xu
				|| yu <= bbox.yl || yl >= bbox.yu ? false : true;
		}
	};

	struct BoundingBoxListItem
	{
		BoundingBox bbox;
		BoundingBoxListItem *next;
	};

	struct Grid
	{
		size_t x_id, y_id;
		BoundingBoxListItem *bbox_list;
		inline bool detect_contact(BoundingBox &bbox)
		{
			for (BoundingBoxListItem *pbox = bbox_list; pbox; pbox = pbox->next)
			{
				if (pbox->bbox.detect_contact(bbox))
					return true;
			}
			return false;
		}

	};

	MemoryUtilities::ItemBuffer<BoundingBoxListItem> bbox_mem;
	inline void add_bbox_to_grid(BoundingBox &bbox, Grid &grid)
	{
		BoundingBoxListItem &bbox_item = *bbox_mem.alloc();
		bbox_item.bbox.copy(bbox);
		bbox_item.next = grid.bbox_list;
		grid.bbox_list = &bbox_item;
	}

	double x0, y0, xn, yn, dx, dy;
	inline bool detect_contact(BoundingBox &bbox)
	{
		return xn <= bbox.xl || x0 >= bbox.xu
			|| yn <= bbox.yl || y0 >= bbox.yu ? false : true;
	}

	Grid *grids;
	size_t x_num, y_num;
	inline Grid *pgrid_by_id(size_t x_id, size_t y_id) { return grids + y_id*x_num + x_id; }

public:
	ParticleContactDetectionGrid() : 
		grids(nullptr), x_num(0), y_num(0) {}
	~ParticleContactDetectionGrid() { clear_grid(); }

	void init_grid(double xl, double xu, double yl, double yu, double hx, double hy);
	void clear_grid(void);

	void reset_grid(void);
	
	size_t add_particle(double x, double y, double vol);
	bool detect_contact(double x, double y, double vol);
};

#endif