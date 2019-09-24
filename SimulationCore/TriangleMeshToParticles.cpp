#include "SimulationCore_pcp.h"

#include "TriangleMeshToParticles.h"

const typename TriangleMeshToParticles::GeneratorFunc
	TriangleMeshToParticles::generator_funcs[] = {
	&FirstOrderGaussPointGenerator,
	&SecondOrderGaussPointGenerator
};

const size_t TriangleMeshToParticles::generator_pcl_num[] = {
	1,
	3
};

const size_t TriangleMeshToParticles::generator_num =
	sizeof(TriangleMeshToParticles::generator_funcs) /
	sizeof(TriangleMeshToParticles::generator_funcs[0]);

void TriangleMeshToParticles::FirstOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, double vol)
{
	Particle &pcl = *add_pcl();
	pcl.x = (p1.x + p2.x + p3.x) / 3.0;
	pcl.y = (p1.y + p2.y + p3.y) / 3.0;
	pcl.vol = vol;
}

void TriangleMeshToParticles::SecondOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, double vol)
{
	vol /= 3.0;
	Particle &pcl1 = *add_pcl();
	pcl1.x = 2.0 / 3.0 * p1.x + 1.0 / 6.0 * p2.x + 1.0 / 6.0 * p3.x;
	pcl1.y = 2.0 / 3.0 * p1.y + 1.0 / 6.0 * p2.y + 1.0 / 6.0 * p3.y;
	pcl1.vol = vol;
	Particle &pcl2 = *add_pcl();
	pcl2.x = 1.0 / 6.0 * p1.x + 2.0 / 3.0 * p2.x + 1.0 / 6.0 * p3.x;
	pcl2.y = 1.0 / 6.0 * p1.y + 2.0 / 3.0 * p2.y + 1.0 / 6.0 * p3.y;
	pcl2.vol = vol;
	Particle &pcl3 = *add_pcl();
	pcl3.x = 1.0 / 6.0 * p1.x + 1.0 / 6.0 * p2.x + 2.0 / 3.0 * p3.x;
	pcl3.y = 1.0 / 6.0 * p1.y + 1.0 / 6.0 * p2.y + 2.0 / 3.0 * p3.y;
	pcl3.vol = vol;
}


void TriangleMeshToParticles::generate_pcls(double max_pcl_area)
{
	size_t elem_num = mesh.get_elem_num();
	TriangleMesh::Element *elems = mesh.get_elems();
	TriangleMesh::Node *nodes = mesh.get_nodes();
	Point elem_n1, elem_n2, elem_n3;
	double elem_area;

	if (max_pcl_area == 0.0)
	{
		particle_buffer.set_page_size(elem_num);
		for (size_t elem_id = 0; elem_id < elem_num; ++elem_id)
		{
			TriangleMesh::Element &elem = elems[elem_id];
			TriangleMesh::Node &n1 = nodes[elem.n1];
			elem_n1.x = n1.x;
			elem_n1.y = n1.y;
			TriangleMesh::Node &n2 = nodes[elem.n2];
			elem_n2.x = n2.x;
			elem_n2.y = n2.y;
			TriangleMesh::Node &n3 = nodes[elem.n3];
			elem_n3.x = n3.x;
			elem_n3.y = n3.y;
			elem_area = abs(((n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (n2.y - n1.y)) / 2.0);
			(this->*cur_generator_func)(elem_n1, elem_n2, elem_n3, elem_area);
		}
	}
	else
	{
		particle_buffer.set_page_size(mesh.get_area() / max_pcl_area / 2.0);
		MemoryUtilities::ItemArray<Point, 2, 22> pt_mem;
		Point *low_pts, *upp_pts;
		double dx21, dy21, dx32, dy32;
		for (size_t elem_id = 0; elem_id < elem_num; ++elem_id)
		{
			TriangleMesh::Element elem = elems[elem_id];
			TriangleMesh::Node &n1 = nodes[elem.n1];
			elem_n1.x = n1.x;
			elem_n1.y = n1.y;
			TriangleMesh::Node &n2 = nodes[elem.n2];
			elem_n2.x = n2.x;
			elem_n2.y = n2.y;
			TriangleMesh::Node &n3 = nodes[elem.n3];
			elem_n3.x = n3.x;
			elem_n3.y = n3.y;
			elem_area = abs(((n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (n2.y - n1.y)) / 2.0);
			size_t div_num = size_t(ceil(sqrt(elem_area / (max_pcl_area * get_pcl_num_per_elem()))));
			if (div_num == 0) div_num = 1;
			elem_area /= double(div_num * div_num);
			dx21 = (elem_n2.x - elem_n1.x) / double(div_num);
			dy21 = (elem_n2.y - elem_n1.y) / double(div_num);
			dx32 = (elem_n3.x - elem_n2.x) / double(div_num);
			dy32 = (elem_n3.y - elem_n2.y) / double(div_num);
			pt_mem.reserve(div_num + div_num + 2);
			upp_pts = pt_mem.get_mem();
			low_pts = upp_pts + div_num + 1;
			low_pts[0].x = elem_n1.x;
			low_pts[0].y = elem_n1.y;
			for (size_t line_id = 0; line_id < div_num; ++line_id)
			{
				Point *tmp = upp_pts;
				upp_pts = low_pts;
				low_pts = tmp;
				low_pts[0].x = upp_pts[0].x + dx21;
				low_pts[0].y = upp_pts[0].y + dy21;
				size_t col_num = line_id + 1;
				for (size_t col_id = 0; col_id < col_num; ++col_id)
				{
					low_pts[col_id + 1].x = low_pts[col_id].x + dx32;
					low_pts[col_id + 1].y = low_pts[col_id].y + dy32;
					(this->*cur_generator_func)(upp_pts[col_id], low_pts[col_id], low_pts[col_id + 1], elem_area);
				}
				for (size_t col_id = 0; col_id < line_id; ++col_id)
				{
					(this->*cur_generator_func)(upp_pts[col_id], low_pts[col_id], upp_pts[col_id + 1], elem_area);
				}
			}
		}
	}
}
