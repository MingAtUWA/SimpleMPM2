#include "TestsWithGL_pcp.h"

#include <ctime>

#include "Model_T2D_CHM_s.h"
#include "DisplayModel_T2D.h"

#include "test_sim_core.h"

extern bool test_AABB_triangle_intersection(double xl, double xu, double yl, double yu,
	double x0, double y0, double x1, double y1, double x2, double y2);

namespace
{
	void print_AABB_tri_intersect(double xl, double xu, double yl, double yu,
		double x0, double y0, double x1, double y1, double x2, double y2)
	{
		printf("AABB (%lf, %lf, %lf, %lf) and triangle (%lf, %lf, %lf, %lf, %lf, %lf) is ",
				xl, xu, yl, yu, x0, y0, x1, y1, x2, y2);
		if (test_AABB_triangle_intersection(xl, xu, yl, yu, x0, y0, x1, y1, x2, y2))
			std::cout << "intersected\n";
		else
			std::cout << "not intersected\n";
	}

	typedef Model_T2D_CHM_s::Element Element;
	typedef Model_T2D_CHM_s::Particle Particle;
}

void test_triangle_searching(void)
{
	//print_AABB_tri_intersect(1.0, 2.0, 1.0, 2.0, 1.51, 0.5, 2.0, 0.0, 2.5, 1.5);
	
	//TriangleMesh tri_mesh;
	//tri_mesh.load_mesh("..\\..\\Asset\\rect.mesh_data");
	////std::cout << "node num: " << tri_mesh.get_node_num() << "\n"
	////		  << "elem num: " << tri_mesh.get_elem_num() << "\n";

	//TriangleMeshToParticles mh_2_pcl(tri_mesh);
	//mh_2_pcl.set_even_div_num(2);
	//mh_2_pcl.set_generator(TriangleMeshToParticles::GeneratorType::RandomlyDistributedPoint);
	//mh_2_pcl.generate_pcls();
	////std::cout << "pcl num: " << mh_2_pcl.get_pcl_num() << "\n";

	//Model_T2D_CHM_s model;
	//// init mesh
	//model.init_mesh(tri_mesh);
	//tri_mesh.clear();
	//// init pcls
	//model.init_pcls(mh_2_pcl, 0.2, 2.0, 1.0, 1000.0, 0.2, 50000.0, 1.0e-4, 1.0);
	//mh_2_pcl.clear();
	//// init bg mesh
	//model.init_bg_mesh(0.05, 0.05);

	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\rect10by15_hole.mesh_data");
	//std::cout << "node num: " << tri_mesh.get_node_num() << "\n"
	//		  << "elem num: " << tri_mesh.get_elem_num() << "\n";

	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	mh_2_pcl.set_generator(TriangleMeshToParticles::GeneratorType::SecondOrderGaussPoint);
	mh_2_pcl.generate_pcls();
	mh_2_pcl.replace_with_grid_points(7.5, 22.5, 7.5, 32.5, 1.0, 1.0);
	//std::cout << "pcl num: " << mh_2_pcl.get_pcl_num() << "\n";

	Model_T2D_CHM_s model;
	// init mesh
	model.init_mesh(tri_mesh);
	tri_mesh.clear();
	// init pcls
	model.init_pcls(mh_2_pcl, 0.2, 2.0, 1.0, 1000.0, 0.2, 50000.0, 1.0e-4, 1.0);
	mh_2_pcl.clear();
	// init bg mesh
	model.init_bg_mesh(1.0, 1.0);
	
	//if (model.is_in_triangle(model.elems[1435], model.pcls[2474]))
	//	std::cout << "true\n";
	//else
	//	std::cout << "false\n";

	//for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	//{
	//	Particle &pcl = model.pcls[p_id];
	//	model.init_pcl_cal_var(pcl);
	//	Element *elem = model.find_in_which_element(pcl);
	//	std::cout << "pcl " << pcl.id << " in elem " 
	//			  << pcl.pe->id << ", " << elem->id << "\n";
	//	//if (pcl.pe->id != elem->id)
	//	//	system("pause");
	//}

	DisplayModel_T2D disp_model;
	disp_model.init_win();
	disp_model.init_model(model);
	disp_model.display(-0.5, 30.5, -0.5, 40.5);

	//system("pause");
}
