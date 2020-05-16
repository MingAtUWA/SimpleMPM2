#ifndef __TEST_POST_PROCESSOR_H__
#define __TEST_POST_PROCESSOR_H__

void test_mpm_rigid_animation_square(void);
void test_mpm_rigid_animation_can(void);
void test_mpm_rigid_animation_can_fric(void);
void test_mpm_rigid_animation_slope_fric(void);
void test_mpm_rigid_animation_bar_compression(void);
void test_animation_chm_s_1d_consolidation(void);

void test_animation_me_s_up_1dbar(void);
void test_animation_chm_s_uup_1d_consolidation(void);

void test_animation_t2d_chm_s_1d_consolidation(void);
void test_color_animation_t2d_chm_s_1d_consolidation(void);

void test_animation_t2d_chm_s_t_bar(void);
void test_color_animation_t2d_chm_s_t_bar(void);

void test_animation_t2d_chm_s_t_bar_coarser(void);
void test_color_animation_t2d_chm_s_t_bar_coarser(void);

void test_color_animation_t2d_chm_s_t_bar_above_ground(void);
void test_color_animation_t2d_chm_s_t_bar_above_ground_geostatic(void);

void test_color_animation_t2d_chm_s_t_bar_real(void);

void test_color_graph(void);

void test_color_animation_t2d_me_s_1d_compression(void);
void test_color_animation_t2d_me_s_t_bar_coarser(void);
void test_color_animation_t2d_me_s_t_bar_above_ground(void);
void test_color_animation_t2d_chm_s_1d_wave(void);

void test_color_animation_t2d_chm_s_geostatic_hdf5(void);
void test_color_animation_t2d_chm_s_geostatic_hdf5_mcc(void);
void test_color_animation_t2d_chm_s_restart_from_geostatic_hdf5_mcc(void);

void test_color_animation_t2d_chm_s_t_bar_real_geostatic(void);
void test_color_animation_t2d_chm_s_t_bar_real_restart(void);

void test_postprocess_t2d_fluid(void);

void test_color_animation_t2d_chm_s_t_bar_conference_geo(void);
void test_color_animation_t2d_chm_s_t_bar_conference_restart(void);
void test_color_animation_t2d_chm_s_t_bar_conference_step2(void);

void test_t2d_chm_geostatic_animation();

void display_s2d_mpm_me_s_geostatic();
void display_s2d_mpm_me_s();

#endif