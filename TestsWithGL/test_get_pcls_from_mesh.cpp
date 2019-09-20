#include "TestsWithGL_pcp.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.hpp"

#include "test_sim_core.h"

void test_get_pcls_from_mesh(void)
{
	double node_coords[] = {
		0.0, 0.0,
		1.0, 0.0,
		0.2, 1.0
	};
	size_t elem_indices[] = {
		0, 1, 2
	};
	TriangleMesh mesh;
	mesh.load_mesh(node_coords, 3, elem_indices, 1);

	struct Particle { double x, y; };
	TriangleMeshToParticles<Particle> mesh_to_pcl(mesh,
		TriangleMeshToParticles<Particle>::GeneratorType::SecondOrderGaussPoint);
	
	size_t pcl_num = mesh_to_pcl.get_pcl_num();
	Particle *pcls = new Particle[pcl_num];
	mesh_to_pcl.get_pcls(pcls);

	std::cout << "pcl_num: " << pcl_num << "\n";
	for (size_t i = 0; i < pcl_num; i++)
		std::cout << "pcl " << i << " - (" << pcls[i].x << ", " << pcls[i].y << ")\n";

	delete[] pcls;

	system("pause");
}
