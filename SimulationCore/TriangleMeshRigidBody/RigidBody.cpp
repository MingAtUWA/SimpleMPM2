#include "SimulationCore_pcp.h"

#include "RigidBody.h"

int RigidBody::load_and_init_mesh(const char *filename, double bg_grid_size)
{
	mesh.load_mesh(filename);
	mesh.get_bounding_box();
	mesh.find_edges_at_boundary();
	mesh.init_bg_grid(bg_grid_size);
	return 0;
}

void RigidBody::get_bounding_box(
		double &x1, double &y1, double &x2, double &y2,
		double &x3, double &y3, double &x4, double &y4,
		double expand_size)
{
	Rect &bd_box = mesh.moved_bounding_box;
	Rect ex_bd_box = { bd_box.xl - expand_size, bd_box.xu + expand_size,
					   bd_box.yl - expand_size, bd_box.yu + expand_size };

	// transform into global coordinates
	//init_transformation();

	double sin_theta_xl = sin_theta * ex_bd_box.xl;
	double cos_theta_xl = cos_theta * ex_bd_box.xl;
	double sin_theta_yl = sin_theta * ex_bd_box.yl;
	double cos_theta_yl = cos_theta * ex_bd_box.yl;
	double sin_theta_xu = sin_theta * ex_bd_box.xu;
	double cos_theta_xu = cos_theta * ex_bd_box.xu;
	double sin_theta_yu = sin_theta * ex_bd_box.yu;
	double cos_theta_yu = cos_theta * ex_bd_box.yu;

	x1 = cos_theta_xl - sin_theta_yl + x;
	y1 = sin_theta_xl + cos_theta_yl + y;
	x2 = cos_theta_xu - sin_theta_yl + x;
	y2 = sin_theta_xu + cos_theta_yl + y;
	x3 = cos_theta_xu - sin_theta_yu + x;
	y3 = sin_theta_xu + cos_theta_yu + y;
	x4 = cos_theta_xl - sin_theta_yu + x;
	y4 = sin_theta_xl + cos_theta_yu + y;
}
