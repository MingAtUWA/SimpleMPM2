#include "TestsWithGL_pcp.h"

#include "ItemArray.hpp"

#include "test_sim_core.h"
#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"
#include "Model_T2D_ME_s.h"
#include "Step_T2D_ME_s_Geostatic.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "DisplayModel_T2D.h"

//#include "test_post_processor.h"
//#include "GA_T2D_ME_s_color.h"

namespace
{

void find_pcl(Model_T2D_ME_s &md)
{
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; pcl_id++)
	{
		Model_T2D_ME_s::Particle &pcl = md.pcls[pcl_id];
		if (pcl.y < 0.02)
		{
			std::cout << pcl.id << " ";
		}
	}
	std::cout << "\n";
}

};

void test_t2d_mpm_me_s_geostatic(void)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\rect.mesh_data");

	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	mh_2_pcl.set_even_div_num(2);
	mh_2_pcl.set_generator(TriangleMeshToParticles::GeneratorType::EvenlyDistributedPoint);
	mh_2_pcl.generate_pcls();

	Model_T2D_ME_s model;
	model.init_mesh(tri_mesh);
	tri_mesh.clear();
	
	model.init_pcls(mh_2_pcl, 20.0, 1000.0, 0.0);
	mh_2_pcl.clear();
	
	model.init_bg_mesh(0.05, 0.05);

	//find_pcl(model);
	//system("pause");
	//return;

	size_t vx_bc_n_id[] = { 0, 3, 15, 16, 17, 18, 19, 20, 21, 22, 23,
							1, 2,  5,  6,  7,  8,  9, 10, 11, 12, 13 };
	model.init_vxs(sizeof(vx_bc_n_id) / sizeof(vx_bc_n_id[0]));
	for (size_t v_id = 0; v_id < model.vx_num; ++v_id)
	{
		VelocityBC &vbc = model.vxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	size_t vy_bc_n_id[] = { 0, 1, 4 };
	model.init_vys(sizeof(vy_bc_n_id) / sizeof(vy_bc_n_id[0]));
	for (size_t v_id = 0; v_id < model.vy_num; ++v_id)
	{
		VelocityBC &vbc = model.vys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	size_t tbc_pcl_id[] = { 132, 133, 168, 169 };
	model.init_tys(sizeof(tbc_pcl_id) / sizeof(tbc_pcl_id[0]));
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM &tbc = model.tys[t_id];
		tbc.pcl_id = tbc_pcl_id[t_id];
		tbc.t = 0.05 * -100.0;
	}

	//MemoryUtilities::ItemArray<GLfloat> pt_array;
	//pt_array.reserve(29 * 3);
	//GLfloat pt_coord;
	//for (size_t n_id = 0; n_id < sizeof(vx_bc_n_id)/sizeof(vx_bc_n_id[0]); ++n_id)
	//{
	//	Model_T2D_ME_s::Node &n = model.nodes[vx_bc_n_id[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}
	//for (size_t n_id = 0; n_id < sizeof(vy_bc_n_id) / sizeof(vy_bc_n_id[0]); ++n_id)
	//{
	//	Model_T2D_ME_s::Node &n = model.nodes[vy_bc_n_id[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}
	//for (size_t p_id = 0; p_id < sizeof(tbc_pcl_id) / sizeof(tbc_pcl_id[0]); ++p_id)
	//{
	//	Model_T2D_ME_s::Particle &pcl = model.pcls[tbc_pcl_id[p_id]];
	//	pt_coord = double(pcl.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(pcl.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}

	//// disp only one point
	//MemoryUtilities::ItemArray<GLfloat> pt_array;
	//GLfloat pt_coord;
	//Model_T2D_ME_s::Particle &pt = model.pcls[5];
	////Model_T2D_ME_s::Node &pt = model.nodes[14];
	//pt_array.reserve(3);
	//pt_coord = double(pt.x);
	//pt_array.add(&pt_coord);
	//pt_coord = double(pt.y);
	//pt_array.add(&pt_coord);
	//pt_coord = 0.0f;
	//pt_array.add(&pt_coord);

	//DisplayModel_T2D disp_model;
	//disp_model.init_win();
	//disp_model.init_model(model);
	////disp_model.init_points(pt_array.get_mem(), pt_array.get_num()/3);
	//disp_model.display(-0.05, 0.25, -0.05, 1.05);
	//return;

	TimeHistoryOutput_ConsoleProgressBar out3;

	Step_T2D_ME_s_Geostatic step;
	step.set_model(model);
	step.set_time(1.0);
	step.set_dtime(1.0e-4);
	step.add_time_history(out3);
	step.solve();

	//system("pause");
}
