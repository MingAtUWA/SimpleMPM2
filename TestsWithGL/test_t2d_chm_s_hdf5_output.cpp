#include "TestsWithGL_pcp.h"

#include "ResultFile_hdf5.h"

#include "Model_T2D_CHM_s.h"
#include "ModelDataOutput_T2D_CHM_s.h"

//#include "Step_T2D_CHM_s_SE.h"
#include "Step_T2D_CHM_s_SE_Geostatic.h"

//#include "TimeHistoryOutput_T2D_CHM_s_SE.h"
#include "TimeHistoryOutput_T2D_CHM_s_SE_Geostatic.h"

void test_t2d_chm_s_hdf5_output(void)
{
	TriangleMesh tri_mesh;
	// elements
	tri_mesh.alloc_elements(2);
	TriangleMesh::Element &e0 = tri_mesh.get_elems()[0];
	e0.id = 0;
	e0.n1 = 0;
	e0.n2 = 1;
	e0.n3 = 2;
	TriangleMesh::Element &e1 = tri_mesh.get_elems()[1];
	e1.id = 1;
	e1.n1 = 0;
	e1.n2 = 2;
	e1.n3 = 3;
	// nodes
	tri_mesh.alloc_nodes(4);
	TriangleMesh::Node &n0 = tri_mesh.get_nodes()[0];
	n0.id = 0;
	n0.x = 0.0;
	n0.y = 0.0;
	TriangleMesh::Node &n1 = tri_mesh.get_nodes()[1];
	n1.id = 1;
	n1.x = 1.0;
	n1.y = 0.0;
	TriangleMesh::Node &n2 = tri_mesh.get_nodes()[2];
	n2.id = 2;
	n2.x = 1.0;
	n2.y = 1.0;
	TriangleMesh::Node &n3 = tri_mesh.get_nodes()[3];
	n3.id = 3;
	n3.x = 0.0;
	n3.y = 1.0;

	Model_T2D_CHM_s model;
	model.set_name("test_model");
	model.init_mesh(tri_mesh);
	model.init_pcls(2, 0.5, 0.5, 2.0, 1.0, 1000.0, 0.25, 1.0e5, 1.0e-5, 1.0);

	ResultFile_hdf5 res_file;
	res_file.create("test_file.hdf5");

	ModelDataOutput_T2D_CHM_s md_out("md_out_test");
	md_out.set_output_time(1.0);
	md_out.set_model(model);
	md_out.set_res_file(res_file);
	md_out.output();

	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic th1("th_out");
	th1.set_interval_num(10);
	th1.set_res_file(res_file);

	Step_T2D_CHM_s_SE_Geostatic step;
	step.set_name("step_test");
	step.set_model(model);
	step.set_time(3.0e-5);
	step.set_dtime(1.0e-5);
	step.add_time_history(th1);
	
	step.solve();

	system("pause");
}