#ifndef __TEST_POST_PROCESSOR_H__
#define __TEST_POST_PROCESSOR_H__

void test_mpm_rigid_animation_square();
void test_mpm_rigid_animation_can();
void test_mpm_rigid_animation_can_fric();
void test_mpm_rigid_animation_slope_fric();
void test_mpm_rigid_animation_bar_compression();
void test_animation_chm_s_1d_consolidation();

void test_animation_me_s_up_1dbar();
void test_animation_chm_s_uup_1d_consolidation();

void test_animation_t2d_chm_s_1d_consolidation();
void test_color_animation_t2d_chm_s_1d_consolidation();

void test_animation_t2d_chm_s_t_bar();
void test_color_animation_t2d_chm_s_t_bar();

void test_animation_t2d_chm_s_t_bar_coarser();
void test_color_animation_t2d_chm_s_t_bar_coarser();

void test_color_graph();

void test_color_animation_t2d_me_s_1d_compression();
void test_color_animation_t2d_me_s_t_bar_coarser();
void test_color_animation_t2d_me_s_t_bar_above_ground();
void test_color_animation_t2d_chm_s_1d_wave();

void test_color_animation_t2d_chm_s_geostatic_hdf5();
void test_color_animation_t2d_chm_s_geostatic_hdf5_mcc();
void test_color_animation_t2d_chm_s_restart_from_geostatic_hdf5_mcc();

void test_postprocess_t2d_fluid();

void test_color_animation_t2d_chm_s_t_bar_conference_geo();
void test_color_animation_t2d_chm_s_t_bar_conference_restart();

void test_t2d_chm_geostatic_animation();

void display_s2d_mpm_me_s_geostatic();
void display_s2d_mpm_me_s();

#endif