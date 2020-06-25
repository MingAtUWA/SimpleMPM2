#ifndef __TEST_SIM_CORE_H__
#define __TEST_SIM_CORE_H__

#include "TriangleMesh.h"

void test_solve_functions();

int display_triangle_mesh(TriangleMesh &tri_mesh, bool disp_tri_mesh, bool dis_bl, bool disp_bg_grid,
						  TriangleMesh::Edge *_edge = nullptr, Point *_pt = nullptr);
void test_matrix_coefficient_set();
void test_cal_stiffness_mat();

void test_triangle_mesh_circle();
void test_triangle_mesh_square();

void test_rigid_body_square();

void test_get_pcls_from_mesh();

void test_mpm_rigidbody_circle();
void test_mpm_rigidbody_square();
void test_mpm_rigidbody_cantilever();
void test_mpm_rigidbody_cantilever_fric();
void test_mpm_rigidbody_bar_compression();

void test_slide_down_frictional_slope();
void test_slide_down_frictional_slope2();

void test_mpm_chm_s_1d_consolidation();

void test_fem_me_s_up_1dbar();
void test_fem_chm_s_1d_consolidation();

void test_imp_mpm_me_s_up_1dbar();
void test_imp_mpm_chm_s_uup_1d_consolidation();

void test_init_pcl_gimp();

// search triangle
void test_triangle_searching();

void test_t2d_mpm_square();
void test_t2d_mpm_chm_s_1d_consolidation();

void test_t2d_mpm_me_s_1d_compression();

void test_t2d_mpm_me_s_t_bar_above_ground();

void test_t2d_mpm_me_s_geostatic();
void test_t2d_mpm_chm_s_geostatic();
void test_t2d_chm_geostatic();

void test_t2d_mpm_chm_s_1d_wave();

void test_t2d_chm_s_hdf5_output();

void test_t2d_chm_s_geostatic_hdf5();
void test_t2d_chm_s_restart_from_geostatic_hdf5();

void test_t2d_chm_s_geostatic_hdf5_mcc();
void test_t2d_chm_s_restart_from_geostatic_hdf5_mcc();

void test_mcc_get_Su();

void test_t2d_fluid();

void test_s2d_mpm_me_s_geostatic();
void test_s2d_mpm_me_s();

void test_s2d_area_distribution();
void test_t2d_area_distribution();
void test_t2d_area_distribution2();

// Modified cam clay
void test_ModifiedCamClay();
void test_mcc_get_Su();
void test_triaxial_test_cam_clay();

void test_t2d_mpm_chm_s_t_bar_conference_geo();
void test_t2d_mpm_chm_s_t_bar_conference_restart1();
void test_t2d_mpm_chm_s_t_bar_conference_restart2();

void test_t2d_chm_restart_1d_consolidation();
void test_t2d_chm_restart_1d_consolidation2();

#endif