#ifndef __TEST_SIM_CORE_H__
#define __TEST_SIM_CORE_H__

#include "TriangleMesh.h"

void test_solve_functions(void);

int display_triangle_mesh(TriangleMesh &tri_mesh, bool disp_tri_mesh, bool dis_bl, bool disp_bg_grid,
						  TriangleMesh::Edge *_edge = nullptr, Point *_pt = nullptr);
void test_matrix_coefficient_set(void);
void test_cal_stiffness_mat(void);

void test_triangle_mesh_circle(void);
void test_triangle_mesh_square(void);

void test_rigid_body_square(void);

void test_get_pcls_from_mesh(void);

void test_mpm_rigidbody_circle(void);
void test_mpm_rigidbody_square(void);
void test_mpm_rigidbody_cantilever(void);
void test_mpm_rigidbody_cantilever_fric(void);
void test_mpm_rigidbody_bar_compression(void);

void test_slide_down_frictional_slope(void);
void test_slide_down_frictional_slope2(void);

void test_mpm_chm_s_1d_consolidation(void);

void test_fem_me_s_up_1dbar(void);
void test_fem_chm_s_1d_consolidation(void);

void test_imp_mpm_me_s_up_1dbar(void);
void test_imp_mpm_chm_s_uup_1d_consolidation(void);

void test_init_pcl_gimp(void);

// search triangle
void test_triangle_searching(void);

void test_t2d_mpm_square(void);
void test_t2d_mpm_chm_s_1d_consolidation(void);
void test_t2d_mpm_chm_s_t_bar(void);
void test_t2d_mpm_chm_s_t_bar_coarser(void);
void test_t2d_mpm_chm_s_t_bar_above_ground(void);

void test_t2d_mpm_me_s_1d_compression(void);
void test_t2d_mpm_me_s_t_bar_coarser(void);

void test_t2d_mpm_chm_s_t_bar_above_ground_geostatic(void);
void test_t2d_mpm_chm_s_t_bar_above_ground_restart(void);
void test_t2d_mpm_me_s_t_bar_above_ground(void);

void test_t2d_mpm_me_s_geostatic(void);
void test_t2d_mpm_chm_s_geostatic(void);

void test_t2d_mpm_chm_s_1d_wave(void);

// Modified cam clay
void test_ModifiedCamClay(void);
void test_mcc_get_Su(void);

void test_t2d_chm_s_hdf5_output(void);

void test_t2d_chm_s_geostatic_hdf5(void);
void test_t2d_chm_s_restart_from_geostatic_hdf5(void);

void test_t2d_chm_s_geostatic_hdf5_mcc(void);
void test_t2d_chm_s_restart_from_geostatic_hdf5_mcc(void);

void test_t2d_mpm_chm_s_t_bar_real_geostatic(void);
void test_t2d_mpm_chm_s_t_bar_real_restart(void);

void test_mcc_get_Su(void);

void test_t2d_fluid(void);

#endif