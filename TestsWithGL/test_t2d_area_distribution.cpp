#include "TestsWithGL_pcp.h"

#include "TriangleMeshToParticles.h"
#include "AdjustParticlesWithTriangleMesh.hpp"
#include "Model_T2D_ME_s.h"
#include "DisplayModel_T2D.h"

#include "test_sim_core.h"

void test_t2d_area_distribution()
{
	double node_coords[] = {
		0.0, 0.0,
		1.0, 0.0,
		2.0, 0.0,
		0.0, 1.0,
		1.0, 1.0,
		2.0, 1.0,
		0.0, 2.0,
		1.0, 2.0,
		2.0, 2.0
	};
	size_t elem_indices[] = {
		0, 1, 3,
		3, 1, 4,
		1, 2, 4,
		4, 2, 5,
		3, 4, 6,
		6, 4, 7,
		4, 5, 7,
		7, 5, 8
	};

	TriangleMesh mesh;
	mesh.load_mesh(node_coords, 9, elem_indices, 8);

	Model_T2D_ME_s model;
	model.init_mesh(mesh);
	model.init_bg_mesh(0.5, 0.5);

	TriangleMeshToParticles t2p(mesh);
	TriangleMeshToParticles::Particle pcl;
	//pcl.x = 0.25;
	//pcl.y = 0.25;
	pcl.x = 1.0;
	pcl.y = 1.0;
	pcl.vol = 0.01;
	t2p.add_pcl(pcl);

	AdjustParticlesWithTriangleMesh<Model_T2D_ME_s> ap_mh(t2p);
	ap_mh.distribute_points_area_to_mesh(model, 0.1, 0.1, 10.0, 4);

	model.init_pcls(t2p, 10.0, 100.0, 0.0);

	DisplayModel_T2D disp_model;
	disp_model.init_win();
	disp_model.init_model(model);
	disp_model.display(-0.05, 2.05, -0.05, 2.05);
}

void test_t2d_area_distribution2()
{
	TriangleMesh mesh;
	mesh.load_mesh("..\\..\\Asset\\rect.mesh_data");

	Model_T2D_ME_s model;
	model.init_mesh(mesh);
	model.init_bg_mesh(0.05, 0.05);

	TriangleMeshToParticles t2p(mesh);
	t2p.generate_grid_points(0.0, 0.2, 0.0, 1.0, 0.025, 0.025);

	AdjustParticlesWithTriangleMesh<Model_T2D_ME_s> ap_mh(t2p);
	ap_mh.distribute_points_area_to_mesh(model, 0.025, 0.025, 0.2, 0.2, 5);

	model.init_pcls(t2p, 10.0, 100.0, 0.0);

	model.sum_vol_for_each_elements();

	MemoryUtilities::ItemArray<GLfloat> pt_array;
	pt_array.reserve(29 * 3);
	GLfloat pt_coord;
	//pt_coord = 0.00859375;
	//pt_array.add(&pt_coord);
	//pt_coord = 0.01640625;
	//pt_array.add(&pt_coord);
	//pt_coord = 0.0f;
	//pt_array.add(&pt_coord);

	DisplayModel_T2D disp_model;
	disp_model.init_win();
	disp_model.init_model(model);
	disp_model.init_points(pt_array.get_mem(), pt_array.get_num() / 3);
	disp_model.display(-0.05, 0.25, -0.05, 1.05);
}