#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <cmath>

struct Point { double x, y; };
struct Rect { double xl, xu, yl, yu; };

inline double distance(Point &p1, Point &p2)
{
	double x_diff = p1.x - p2.x;
	double y_diff = p1.y - p2.y;
	return sqrt(x_diff * x_diff + y_diff * y_diff);
}

inline double distance(Rect &rec, Point &p)
{
	double cx, cy, x_diff, y_diff;
	cx = p.x < rec.xu ? p.x : rec.xu;
	cx = rec.xl > cx ? rec.xl : cx;
	cy = p.y < rec.yu ? p.y : rec.yu;
	cy = rec.yl > cy ? rec.yl : cy;
	x_diff = p.x - cx;
	y_diff = p.y - cy;
	return sqrt(x_diff * x_diff + y_diff * y_diff);
}

bool clip_line(double xl, double xu, double yl, double yu,
			   double &x1, double &y1, double &x2, double &y2);

// whether AABB and triangle intersect
// return true:  intersect
// return false: intersect
bool test_AABB_triangle_intersection(double xl, double xu, double yl, double yu,
	double x0, double y0, double x1, double y1, double x2, double y2);

#endif