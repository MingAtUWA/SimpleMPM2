#include "TestsWithGL_pcp.h"

#include "Model_S2D_CHM_s.h"

#include "test_sim_core.h"

void disp_pcl_vars(Model_S2D_CHM_s::Particle &pcl)
{
	std::cout << "Particle " << pcl.index << ": (x, y, vol)\n";
	std::cout << "  s - " << pcl.var.x << ", " << pcl.var.y
			  << ", " << pcl.var.vol << "\n";
	for (size_t i = 0; i < pcl.elem_num; i++)
	{
		Model_S2D_CHM_s::ParticleCalVar &pcl_var = pcl.vars[i];
		std::cout << "  " << i << " - " << pcl_var.x << ", "
				  << pcl_var.y << ", " << pcl_var.vol << "\n";
	}
}

static Model_S2D_CHM_s model;

void test_pcl(double x, double y, double vol)
{
	static size_t cur_id = 0;
	Model_S2D_CHM_s::Particle &pcl = model.pcls[cur_id++];
	pcl.x = x;
	pcl.y = y;
	pcl.m_s = vol;
	model.init_pcl_GIMP(pcl);
	disp_pcl_vars(pcl);
}

void test_init_pcl_gimp(void)
{
	model.init_mesh(0.5, 10, 10, 0.0, 0.0, 0.05);
	model.init_pcl(30, 0.0, 1.0, 1.0, 1.0, 100.0, 0.0, 5000.0, 1.0e-3, 1.0);

	test_pcl(0.25, 0.25, 0.2);

	// 2.3, 2.9 - 1.7, 2.3
	test_pcl(2.6, 2.0, 0.36);

	// 2.3, 3.0 - 1.7, 2.4
	test_pcl(2.65, 2.05, 0.49);

	// 2.5, 3.0 - 2.0, 2.5
	test_pcl(2.75, 2.25, 0.25);

	// 2.5, 3.0 - 1.75, 2.25
	test_pcl(2.75, 2.0, 0.25);

	// 2.5, 3.0 - 1.75, 2.25
	test_pcl(2.75, 2.0, 0.16);

	// 2.3, 2.7 - 2.05, 2.45
	test_pcl(2.6, 2.25, 0.16);

	// 2.3, 2.7 - 2.5, 2.9
	test_pcl(2.4, 2.7, 0.16);

	// 2.2, 2.6 - 2.4, 2.8
	test_pcl(2.4, 2.6, 0.16);

	// 2.4, 3.1 - 1.9, 2.6
	test_pcl(2.75, 2.25, 0.49);

	// 2.0, 3.5 - 1.5, 3.0
	test_pcl(2.8, 2.3, 2.56);

	// 2.0, 3.5 - 1.5, 3.0
	test_pcl(2.79, 2.29, 2.56);

	// at boundary 2.5, 3.0 - 0.0, 0.3
	test_pcl(2.75, 0.1, 0.16);

	// at boundary 2.5, 3.0 - 0.0, 0.3
	test_pcl(4.9, 0.1, 0.16);

	// at boundary 2.5, 3.0 - 0.0, 0.3
	test_pcl(3.75, 0.1, 2.56);

	// at boundary 2.5, 3.0 - 0.0, 0.3
	test_pcl(3.75, 0.799999, 2.56);

	system("pause");
}
