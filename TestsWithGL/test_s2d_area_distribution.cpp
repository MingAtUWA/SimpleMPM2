#include "TestsWithGL_pcp.h"

#include "Step_S2D_ME_s_Geostatic.h"

#include "test_sim_core.h"

namespace
{
	void print_area(Model_S2D_ME_s &md)
	{
		Model_S2D_ME_s::Element *elems = md.get_elems();
		for (size_t i = md.get_elem_y_num(); i > 0; --i)
		{
			for (size_t j = 0; j < md.get_elem_x_num(); ++j)
			{
				std::cout << elems[(i-1) * md.get_elem_x_num() + j].ve_vol << " ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}

	void reset_elem_area(Model_S2D_ME_s &md)
	{
		for (size_t i = 0; i < md.get_elem_num(); i++)
		{
			Model_S2D_ME_s::Element &e = md.get_elems()[i];
			e.ve_vol = 0.0;
		}
	}
}

void test_s2d_area_distribution()
{
	typedef Model_S2D_ME_s::Element Element;
	typedef Model_S2D_ME_s::Particle Particle;
	
	Model_S2D_ME_s md;

	size_t elem_x_num = 4;
	size_t elem_y_num = 4;
	md.init_mesh(0.0, 0.0, 2.0,2.0, elem_x_num, elem_y_num);

	md.alloc_pcls(10);
	Particle *pcls = md.get_pcls();

	reset_elem_area(md);
	Particle &pcl0 = pcls[0];
	pcl0.x = 1.25;
	pcl0.y = 1.25;
	pcl0.m = 0.09;
	pcl0.density = 1.0;
	pcl0.ar = 1.0;
	md.distribute_pcl_mat_to_elems(pcl0);
	print_area(md);

	reset_elem_area(md);
	Particle &pcl1 = pcls[1];
	pcl1.x = 1.5;
	pcl1.y = 1.25;
	pcl1.m = 0.15;
	pcl1.density = 1.0;
	pcl1.ar = 2.0;
	md.distribute_pcl_mat_to_elems(pcl1);
	print_area(md);

	reset_elem_area(md);
	Particle &pcl2 = pcls[2];
	pcl2.x = 0.75;
	pcl2.y = 1.5;
	pcl2.m = 0.15;
	pcl2.density = 1.0;
	pcl2.ar = 0.5;
	md.distribute_pcl_mat_to_elems(pcl2);
	print_area(md);

	reset_elem_area(md);
	Particle &pcl3 = pcls[3];
	pcl3.x = 0.5;
	pcl3.y = 1.0;
	pcl3.m = 0.16;
	pcl3.density = 1.0;
	pcl3.ar = 1.0;
	md.distribute_pcl_mat_to_elems(pcl3);
	print_area(md);

	reset_elem_area(md);
	Particle &pcl4 = pcls[4];
	pcl4.x = 0.75;
	pcl4.y = 0.55;
	pcl4.m = 0.49;
	pcl4.density = 1.0;
	pcl4.ar = 1.0;
	md.distribute_pcl_mat_to_elems(pcl4);
	print_area(md);

	reset_elem_area(md);
	Particle &pcl5 = pcls[5];
	pcl5.x = 1.75;
	pcl5.y = 0.75;
	pcl5.m = 0.3;
	pcl5.density = 1.0;
	pcl5.ar = 0.3;
	md.distribute_pcl_mat_to_elems(pcl5);
	print_area(md);

	reset_elem_area(md);
	Particle &pcl6 = pcls[6];
	pcl6.x = 1.0;
	pcl6.y = 1.25;
	pcl6.m = 0.45;
	pcl6.density = 1.0;
	pcl6.ar = 5.0;
	md.distribute_pcl_mat_to_elems(pcl6);
	print_area(md);
}
