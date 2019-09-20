#ifndef __TEST_SIM_CORE_H__
#define __TEST_SIM_CORE_H__

#include "TriangleMesh.h"

int display_triangle_mesh(TriangleMesh &tri_mesh, bool disp_tri_mesh, bool dis_bl, bool disp_bg_grid,
						  TriangleMesh::Edge *_edge = nullptr, Point *_pt = nullptr);

void test_triangle_mesh_circle(void);
void test_triangle_mesh_square(void);

void test_rigid_body_square(void);

void test_mpm_rigidbody_circle(void);
void test_mpm_rigidbody_square(void);
void test_mpm_rigidbody_cantilever(void);
void test_mpm_rigidbody_cantilever_fric(void);

void test_get_pcls_from_mesh(void);

#endif