#include "TestsWithGL_pcp.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"

#include "test_sim_core.h"

void test_get_pcls_from_mesh(void)
{
	double node_coords[] = {
		0.0, 0.0,
		1.0, 0.0,
		0.0, 1.0,
		1.0, 1.0
	};
	size_t elem_indices[] = {
		0, 1, 2,
		2, 1, 3
	};
	TriangleMesh mesh;
	mesh.load_mesh(node_coords, 4, elem_indices, 2);

	struct Particle { double x, y; };
	TriangleMeshToParticles mesh_to_pcl(mesh,
		//TriangleMeshToParticles::GeneratorType::FirstOrderGuassPoint);
		TriangleMeshToParticles::GeneratorType::SecondOrderGaussPoint);

	mesh_to_pcl.generate_pcls(0.15);
	size_t pcl_num = mesh_to_pcl.get_pcl_num();

	std::cout << "pcl_num: " << pcl_num << "\n";
	size_t pcl_id = 0;
	for (TriangleMeshToParticles::Particle *ppcl = mesh_to_pcl.first();
		ppcl; ppcl = mesh_to_pcl.next(ppcl))
	{
		std::cout << "pcl " << pcl_id
				  << " - (" << ppcl->x << ", " << ppcl->y << ") "
				  << ppcl->vol << "\n";
		++pcl_id;
	}
	
	system("pause");
}
