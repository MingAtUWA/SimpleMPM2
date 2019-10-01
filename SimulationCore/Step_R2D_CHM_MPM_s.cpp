#include "SimulationCore_pcp.h"

#include "Step_R2D_CHM_MPM_s.h"

Step_R2D_CHM_MPM_s::Step_R2D_CHM_MPM_s() :
	Step(&solve_substep_R2D_CHM_MPM_s),
	model(nullptr) {}

Step_R2D_CHM_MPM_s::~Step_R2D_CHM_MPM_s() {}

int Step_R2D_CHM_MPM_s::init()
{
	Particle_R2D_CHM_s *ppcl;

	if (is_first_step)
	{
		for (size_t i = 0; i < model->pcl_num; i++)
		{
			ppcl = model->pcls + i;
			ppcl->is_in_mesh = true;
			ppcl->elem = nullptr;
		}
	}

	for (size_t i = 0; i < model->pcl_num; i++)
	{
		ppcl = model->pcls + i;
		ppcl->x_ori = ppcl->x;
		ppcl->y_ori = ppcl->y;
		ppcl->ux_s = 0.0;
		ppcl->uy_s = 0.0;
		ppcl->ux_f = 0.0;
		ppcl->uy_f = 0.0;
	}

	return 0;
}

int Step_R2D_CHM_MPM_s::finalize() { return 0; }

int solve_substep_R2D_CHM_MPM_s(void *_self)
{
	Step_R2D_CHM_MPM_s *self = (Step_R2D_CHM_MPM_s *)_self;
	Model_R2D_CHM_MPM_s *model = self->model;
	Particle_R2D_CHM_s *ppcl;
	Element_R2D_CHM_MPM_s *pelem;
	Node_R2D_CHM_s *pn, *pn1, *pn2, *pn3, *pn4;

	// init nodes
	for (size_t i = 0; i < model->node_num; i++)
	{
		pn = model->nodes + i;
		pn->cal_flag = 0;
		
		// for mapping from nodes
		pn->m_s = 0.0;
		pn->mmx_s = 0.0;
		pn->mmy_s = 0.0;
		pn->fx_kin_f = 0.0;
		pn->fy_kin_f = 0.0;
		pn->fx_ext_m = 0.0;
		pn->fy_ext_m = 0.0;
		pn->fx_int_m = 0.0;
		pn->fy_int_m = 0.0;
		// for mapping back to particles
		pn->ax_s = 0.0;
		pn->ay_s = 0.0;
		pn->vx_s = 0.0;
		pn->vy_s = 0.0;
		pn->dux_s = 0.0;
		pn->duy_s = 0.0;

		// for mapping from nodes
		pn->m_tf = 0.0;
		pn->mmx_tf = 0.0;
		pn->mmy_tf = 0.0;
		pn->fx_ext_tf = 0.0;
		pn->fy_ext_tf = 0.0;
		pn->fx_int_tf = 0.0;
		pn->fy_int_tf = 0.0;
		pn->fx_drag_tf = 0.0;
		pn->fy_drag_tf = 0.0;
		// for mapping back to particles
		pn->ax_f = 0.0;
		pn->ay_f = 0.0;
		pn->vx_f = 0.0;
		pn->vy_f = 0.0;
		pn->dux_f = 0.0;
		pn->duy_f = 0.0;
	}

	// init particles
	for (size_t i = 0; i < model->pcl_num; i++)
	{
		ppcl = model->pcls + i;
		if (ppcl->is_in_mesh)
		{
			pelem = model->find_in_which_element(ppcl->x, ppcl->y, ppcl->elem);
			ppcl->elem = pelem;
			if (pelem)
			{
				// init variables on particles
				ppcl->k_div_miu = ppcl->k / ppcl->miu;
				
				// cal shape functions and their derivatives
				model->cal_shape_function(ppcl);
				
				// init the four nodes associated with this particle
				Get_Nodes_Of_Element_R2D(pelem, model, pn1, pn2, pn3, pn4);
				ppcl->node1 = pn1;
				pn1->cal_flag = 1;
				ppcl->node2 = pn2;
				pn2->cal_flag = 1;
				ppcl->node3 = pn3;
				pn3->cal_flag = 1;
				ppcl->node4 = pn4;
				pn4->cal_flag = 1;
			}
			else
			{
				ppcl->is_in_mesh = false;
			}
		}
	}

	// map variables to node and cal internal force
	for (size_t i = 0; i < model->pcl_num; i++)
	{
		ppcl = model->pcls + i;
		if (ppcl->is_in_mesh)
		{
			// ------------------- node 1 -------------------
			pn1 = ppcl->node1;
			// mixture phase
			pn1->m_s   += ppcl->N1 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vol;
			pn1->mmx_s += ppcl->N1 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vx_s * ppcl->vol;
			pn1->mmy_s += ppcl->N1 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vy_s * ppcl->vol;
			pn1->fx_int_m += (ppcl->dN1_dx * (ppcl->s11 - ppcl->p) + ppcl->dN1_dy * ppcl->s12) * ppcl->vol;
			pn1->fy_int_m += (ppcl->dN1_dx * ppcl->s12 + ppcl->dN1_dy * (ppcl->s22 - ppcl->p)) * ppcl->vol;
			// fluid phase
			pn1->m_tf   += ppcl->N1 * ppcl->density_f * ppcl->vol;
			pn1->mmx_tf += ppcl->N1 * ppcl->density_f * ppcl->vx_f * ppcl->vol;
			pn1->mmy_tf += ppcl->N1 * ppcl->density_f * ppcl->vy_f * ppcl->vol;
			pn1->fx_int_tf  += (ppcl->dN1_dx * -ppcl->p) * ppcl->vol;
			pn1->fy_int_tf  += (ppcl->dN1_dy * -ppcl->p) * ppcl->vol;
			pn1->fx_drag_tf += ppcl->N1 * ppcl->n / ppcl->k_div_miu * (ppcl->vx_f - ppcl->vx_s) * ppcl->vol;
			pn1->fy_drag_tf += ppcl->N1 * ppcl->n / ppcl->k_div_miu * (ppcl->vy_f - ppcl->vy_s) * ppcl->vol;

			// ------------------- node 2 -------------------
			pn2 = ppcl->node2;
			// mixture phase
			pn2->m_s   += ppcl->N2 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vol;
			pn2->mmx_s += ppcl->N2 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vx_s * ppcl->vol;
			pn2->mmy_s += ppcl->N2 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vy_s * ppcl->vol;
			pn2->fx_int_m += (ppcl->dN2_dx * (ppcl->s11 - ppcl->p) + ppcl->dN2_dy * ppcl->s12) * ppcl->vol;
			pn2->fy_int_m += (ppcl->dN2_dx * ppcl->s12 + ppcl->dN2_dy * (ppcl->s22 - ppcl->p)) * ppcl->vol;
			// fluid phase
			pn2->m_tf   += ppcl->N2 * ppcl->density_f * ppcl->vol;
			pn2->mmx_tf += ppcl->N2 * ppcl->density_f * ppcl->vx_f * ppcl->vol;
			pn2->mmy_tf += ppcl->N2 * ppcl->density_f * ppcl->vy_f * ppcl->vol;
			pn2->fx_int_tf += (ppcl->dN2_dx * -ppcl->p) * ppcl->vol;
			pn2->fy_int_tf += (ppcl->dN2_dy * -ppcl->p) * ppcl->vol;
			pn2->fx_drag_tf += ppcl->N2 * ppcl->n / ppcl->k_div_miu * (ppcl->vx_f - ppcl->vx_s) * ppcl->vol;
			pn2->fy_drag_tf += ppcl->N2 * ppcl->n / ppcl->k_div_miu * (ppcl->vy_f - ppcl->vy_s) * ppcl->vol;
			
			// ------------------- node 3 -------------------
			pn3 = ppcl->node3;
			// mixture phase
			pn3->m_s   += ppcl->N3 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vol;
			pn3->mmx_s += ppcl->N3 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vx_s * ppcl->vol;
			pn3->mmy_s += ppcl->N3 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vy_s * ppcl->vol;
			pn3->fx_int_m += (ppcl->dN3_dx * (ppcl->s11 - ppcl->p) + ppcl->dN3_dy * ppcl->s12) * ppcl->vol;
			pn3->fy_int_m += (ppcl->dN3_dx * ppcl->s12 + ppcl->dN3_dy * (ppcl->s22 - ppcl->p)) * ppcl->vol;
			// fluid phase
			pn3->m_tf   += ppcl->N3 * ppcl->density_f * ppcl->vol;
			pn3->mmx_tf += ppcl->N3 * ppcl->density_f * ppcl->vx_f * ppcl->vol;
			pn3->mmy_tf += ppcl->N3 * ppcl->density_f * ppcl->vy_f * ppcl->vol;
			pn3->fx_int_tf += (ppcl->dN3_dx * -ppcl->p) * ppcl->vol;
			pn3->fy_int_tf += (ppcl->dN3_dy * -ppcl->p) * ppcl->vol;
			pn3->fx_drag_tf += ppcl->N3 * ppcl->n / ppcl->k_div_miu * (ppcl->vx_f - ppcl->vx_s) * ppcl->vol;
			pn3->fy_drag_tf += ppcl->N3 * ppcl->n / ppcl->k_div_miu * (ppcl->vy_f - ppcl->vy_s) * ppcl->vol;

			// ------------------- node 4 -------------------
			pn4 = ppcl->node4;
			// mixture phase
			pn4->m_s   += ppcl->N4 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vol;
			pn4->mmx_s += ppcl->N4 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vx_s * ppcl->vol;
			pn4->mmy_s += ppcl->N4 * (1.0 - ppcl->n) * ppcl->density_s * ppcl->vy_s * ppcl->vol;
			pn4->fx_int_m += (ppcl->dN4_dx * (ppcl->s11 - ppcl->p) + ppcl->dN4_dy * ppcl->s12) * ppcl->vol;
			pn4->fy_int_m += (ppcl->dN4_dx * ppcl->s12 + ppcl->dN4_dy * (ppcl->s22 - ppcl->p)) * ppcl->vol;
			// fluid phase
			pn4->m_tf   += ppcl->N4 * ppcl->density_f * ppcl->vol;
			pn4->mmx_tf += ppcl->N4 * ppcl->density_f * ppcl->vx_f * ppcl->vol;
			pn4->mmy_tf += ppcl->N4 * ppcl->density_f * ppcl->vy_f * ppcl->vol;
			pn4->fx_int_tf += (ppcl->dN4_dx * -ppcl->p) * ppcl->vol;
			pn4->fy_int_tf += (ppcl->dN4_dy * -ppcl->p) * ppcl->vol;
			pn4->fx_drag_tf += ppcl->N4 * ppcl->n / ppcl->k_div_miu * (ppcl->vx_f - ppcl->vx_s) * ppcl->vol;
			pn4->fy_drag_tf += ppcl->N4 * ppcl->n / ppcl->k_div_miu * (ppcl->vy_f - ppcl->vy_s) * ppcl->vol;
		}
	}

	// body force
	double bf_m, bf_tf;
	for (size_t i = 0; i < model->bfx_num; i++)
	{
		ppcl  = model->pcls + model->bfxs[i].pcl_id;
		// body force on particle
		bf_m  = ppcl->vol * ((1.0 - ppcl->n) * ppcl->density_s + ppcl->n * ppcl->density_f) * model->bfxs[i].bf;
		bf_tf = ppcl->vol * ppcl->density_f * model->bfxs[i].bf;
		// node 1
		pn1 = ppcl->node1;
		pn1->fx_ext_m  += ppcl->N1 * bf_m;
		pn1->fx_ext_tf += ppcl->N1 * bf_tf;
		// node 2
		pn2 = ppcl->node2;
		pn2->fx_ext_m  += ppcl->N2 * bf_m;
		pn2->fx_ext_tf += ppcl->N2 * bf_tf;
		// node 3
		pn3 = ppcl->node3;
		pn3->fx_ext_m  += ppcl->N3 * bf_m;
		pn3->fx_ext_tf += ppcl->N3 * bf_tf;
		// node 4
		pn4 = ppcl->node4;
		pn4->fx_ext_m  += ppcl->N4 * bf_m;
		pn4->fx_ext_tf += ppcl->N4 * bf_tf;
	}
	for (size_t i = 0; i < model->bfy_num; i++)
	{
		ppcl  = model->pcls + model->bfys[i].pcl_id;
		// body force on particle
		bf_m  = ppcl->vol * ((1.0 - ppcl->n) * ppcl->density_s + ppcl->n * ppcl->density_f) * model->bfys[i].bf;
		bf_tf = ppcl->vol * ppcl->density_f * model->bfys[i].bf;
		// node 1
		pn1 = ppcl->node1;
		pn1->fy_ext_m  += ppcl->N1 * bf_m;
		pn1->fy_ext_tf += ppcl->N1 * bf_tf;
		// node 2
		pn2 = ppcl->node2;
		pn2->fy_ext_m  += ppcl->N2 * bf_m;
		pn2->fy_ext_tf += ppcl->N2 * bf_tf;
		// node 3
		pn3 = ppcl->node3;
		pn3->fy_ext_m  += ppcl->N3 * bf_m;
		pn3->fy_ext_tf += ppcl->N3 * bf_tf;
		// node 4
		pn4 = ppcl->node4;
		pn4->fy_ext_m  += ppcl->N4 * bf_m;
		pn4->fy_ext_tf += ppcl->N4 * bf_tf;
	}

	// surface force
	for (size_t i = 0; i < model->tx_bc_num; i++)
	{
		ppcl = model->pcls + model->tx_bcs[i].pcl_id;
		// node 1
		pn1 = ppcl->node1;
		pn1->fx_ext_m += ppcl->N1 * model->tx_bcs[i].t;
		// node 2
		pn2 = ppcl->node2;
		pn2->fx_ext_m += ppcl->N2 * model->tx_bcs[i].t;
		// node 3
		pn3 = ppcl->node3;
		pn3->fx_ext_m += ppcl->N3 * model->tx_bcs[i].t;
		// node 4
		pn4 = ppcl->node4;
		pn4->fx_ext_m += ppcl->N4 * model->tx_bcs[i].t;
	}
	for (size_t i = 0; i < model->ty_bc_num; i++)
	{
		ppcl = model->pcls + model->ty_bcs[i].pcl_id;
		// node 1
		pn1 = ppcl->node1;
		pn1->fy_ext_m += ppcl->N1 * model->ty_bcs[i].t;
		// node 2
		pn2 = ppcl->node2;
		pn2->fy_ext_m += ppcl->N2 * model->ty_bcs[i].t;
		// node 3
		pn3 = ppcl->node3;
		pn3->fy_ext_m += ppcl->N3 * model->ty_bcs[i].t;
		// node 4
		pn4 = ppcl->node4;
		pn4->fy_ext_m += ppcl->N4 * model->ty_bcs[i].t;
	}
	// pore pressure force...

	// update nodal acceleration of fluid pahse
	for (size_t i = 0; i < model->node_num; i++)
	{
		pn = model->nodes + i;
		if (pn->cal_flag)
		{
			pn->ax_f = (pn->fx_ext_tf - pn->fx_int_tf - pn->fx_drag_tf) / pn->m_tf;
			pn->ay_f = (pn->fy_ext_tf - pn->fy_int_tf - pn->fy_drag_tf) / pn->m_tf;
		}
		//std::cout << "2: fint_tf " << pn->fx_int_tf << "\n";
	}
	for (size_t i = 0; i < model->ax_f_bc_num; i++)
	{
		pn = model->nodes + model->ax_f_bcs[i].node_id;
		pn->ax_f = model->ax_f_bcs[i].a;
	}
	for (size_t i = 0; i < model->ay_f_bc_num; i++)
	{
		pn = model->nodes + model->ay_f_bcs[i].node_id;
		pn->ay_f = model->ay_f_bcs[i].a;
	}

	// update nodal momentum of fluid phase
	for (size_t i = 0; i < model->node_num; i++)
	{
		pn = model->nodes + i;
		if (pn->cal_flag)
		{
			pn->vx_f  = pn->mmx_tf / pn->m_tf;
			pn->vx_f += pn->ax_f * self->dt;
			pn->vy_f  = pn->mmy_tf / pn->m_tf;
			pn->vy_f += pn->ay_f * self->dt;
		}
	}
	// apply velocity boundary conditions of fluid phase
	for (size_t i = 0; i < model->vx_f_bc_num; i++)
	{
		pn = model->nodes + model->vx_f_bcs[i].node_id;
		pn->vx_f = model->vx_f_bcs[i].v;
		pn->ax_f = 0.0;
	}
	for (size_t i = 0; i < model->vy_f_bc_num; i++)
	{
		pn = model->nodes + model->vy_f_bcs[i].node_id;
		pn->vy_f = model->vy_f_bcs[i].v;
		pn->ay_f = 0.0;
	}

	// calculate the inertial term of fluid in mixture formulation
	double pcl_ax_f, pcl_ay_f;
	double pcl_max_f, pcl_may_f;
	for (size_t i = 0; i < model->pcl_num; i++)
	{
		ppcl = model->pcls + i;
		pn1 = ppcl->node1;
		pn2 = ppcl->node2;
		pn3 = ppcl->node3;
		pn4 = ppcl->node4;
		// particle acceleration
		pcl_ax_f = ppcl->N1 * pn1->ax_f + ppcl->N2 * pn2->ax_f 
				 + ppcl->N3 * pn3->ax_f + ppcl->N4 * pn4->ax_f;
		pcl_ay_f = ppcl->N1 * pn1->ay_f + ppcl->N2 * pn2->ay_f
				 + ppcl->N3 * pn3->ay_f + ppcl->N4 * pn4->ay_f;
		
		pcl_max_f = ppcl->n * ppcl->density_f * pcl_ax_f * ppcl->vol;
		pcl_may_f = ppcl->n * ppcl->density_f * pcl_ay_f * ppcl->vol;
		// node 1
		pn1 = ppcl->node1;
		pn1->fx_kin_f += ppcl->N1 * pcl_max_f;
		pn1->fy_kin_f += ppcl->N1 * pcl_may_f;
		// node 2
		pn2 = ppcl->node2;
		pn2->fx_kin_f += ppcl->N2 * pcl_max_f;
		pn2->fy_kin_f += ppcl->N2 * pcl_may_f;
		// node 3
		pn3 = ppcl->node3;
		pn3->fx_kin_f += ppcl->N3 * pcl_max_f;
		pn3->fy_kin_f += ppcl->N3 * pcl_may_f;
		// node 4
		pn4 = ppcl->node4;
		pn4->fx_kin_f += ppcl->N4 * pcl_max_f;
		pn4->fy_kin_f += ppcl->N4 * pcl_may_f;
	}

	// update nodal velocity of solid phase
	for (size_t i = 0; i < model->node_num; i++)
	{
		pn = model->nodes + i;
		if (pn->cal_flag)
		{
			pn->ax_s = (pn->fx_ext_m - pn->fx_int_m - pn->fx_kin_f) / pn->m_s;
			pn->ay_s = (pn->fy_ext_m - pn->fy_int_m - pn->fy_kin_f) / pn->m_s;
		}
	}
	// apply acceleration boundary conditions
	for (size_t i = 0; i < model->ax_s_bc_num; i++)
	{
		pn = model->nodes + model->ax_s_bcs[i].node_id;
		pn->ax_s = model->ax_s_bcs[i].a;
	}
	for (size_t i = 0; i < model->ay_s_bc_num; i++)
	{
		pn = model->nodes + model->ay_s_bcs[i].node_id;
		pn->ay_s = model->ay_s_bcs[i].a;
	}
	
	// update nodal momentum of fluid pahse
	for (size_t i = 0; i < model->node_num; i++)
	{
		pn = model->nodes + i;
		if (pn->cal_flag)
		{
			pn->vx_s  = pn->mmx_s / pn->m_s;
			pn->vx_s += pn->ax_s * self->dt;
			pn->vy_s  = pn->mmy_s / pn->m_s;
			pn->vy_s += pn->ay_s * self->dt;
		}
	}
	// apply velocity boundary conditions of solid phase
	for (size_t i = 0; i < model->vx_s_bc_num; i++)
	{
		pn = model->nodes + model->vx_s_bcs[i].node_id;
		pn->vx_s = model->vx_s_bcs[i].v;
		pn->ax_s = 0.0;
	}
	for (size_t i = 0; i < model->vy_s_bc_num; i++)
	{
		pn = model->nodes + model->vy_s_bcs[i].node_id;
		pn->vy_s = model->vy_s_bcs[i].v;
		pn->ay_s = 0.0;
	}

	// update displacement increment of both phases
	for (size_t i = 0; i < model->node_num; i++)
	{
		pn = model->nodes + i;
		if (pn->cal_flag)
		{
			// solid phase
			pn->dux_s = pn->vx_s * self->dt;
			pn->duy_s = pn->vy_s * self->dt;
			// fluid phase
			pn->dux_f = pn->vx_f * self->dt;
			pn->duy_f = pn->vy_f * self->dt;
			//if (i == 3)
			//	std::cout << "1 vxs: " << pn->vx_s << " vys: " << pn->vy_s
			//			  <<  " vxf: " << pn->vx_f << " vyf: " << pn->vy_f << "\n";
		}
	}

	// map variables back to and update variables particles
	double N1_tmp, N2_tmp, N3_tmp, N4_tmp;
	double E_tmp;
	double de11_s, de22_s, de12_s, dw12;
	double ds11, ds22, ds12;
	double de_vol_s, de_vol_f;
	for (size_t i = 0; i < model->pcl_num; i++)
	{
		ppcl = model->pcls + i;
		if (ppcl->is_in_mesh)
		{
			pn1 = ppcl->node1;
			pn2 = ppcl->node2;
			pn3 = ppcl->node3;
			pn4 = ppcl->node4;
			N1_tmp = ppcl->N1;
			N2_tmp = ppcl->N2;
			N3_tmp = ppcl->N3;
			N4_tmp = ppcl->N4;

			// velocity
			ppcl->vx_s += (pn1->ax_s * N1_tmp + pn2->ax_s * N2_tmp
						 + pn3->ax_s * N3_tmp + pn4->ax_s * N4_tmp) * self->dt;
			ppcl->vy_s += (pn1->ay_s * N1_tmp + pn2->ay_s * N2_tmp
						 + pn3->ay_s * N3_tmp + pn4->ay_s * N4_tmp) * self->dt;
			ppcl->vx_f += (pn1->ax_f * N1_tmp + pn2->ax_f * N2_tmp
						 + pn3->ax_f * N3_tmp + pn4->ax_f * N4_tmp) * self->dt;
			ppcl->vy_f += (pn1->ay_f * N1_tmp + pn2->ay_f * N2_tmp
						 + pn3->ay_f * N3_tmp + pn4->ay_f * N4_tmp) * self->dt;
			

			// displacement
			ppcl->ux_s += pn1->dux_s * N1_tmp + pn2->dux_s * N2_tmp
						+ pn3->dux_s * N3_tmp + pn4->dux_s * N4_tmp;
			ppcl->uy_s += pn1->duy_s * N1_tmp + pn2->duy_s * N2_tmp
						+ pn3->duy_s * N3_tmp + pn4->duy_s * N4_tmp;
			ppcl->ux_f += pn1->dux_f * N1_tmp + pn2->dux_f * N2_tmp
						+ pn3->dux_f * N3_tmp + pn4->dux_f * N4_tmp;
			ppcl->uy_f += pn1->duy_f * N1_tmp + pn2->duy_f * N2_tmp
						+ pn3->duy_f * N3_tmp + pn4->duy_f * N4_tmp;

			// update position
			ppcl->x = ppcl->x_ori + ppcl->ux_s;
			ppcl->y = ppcl->y_ori + ppcl->uy_s;

			// strain increment
			de11_s = pn1->dux_s * ppcl->dN1_dx + pn2->dux_s * ppcl->dN2_dx
				   + pn3->dux_s * ppcl->dN3_dx + pn4->dux_s * ppcl->dN4_dx;
			de22_s = pn1->duy_s * ppcl->dN1_dy + pn2->duy_s * ppcl->dN2_dy
				   + pn3->duy_s * ppcl->dN3_dy + pn4->duy_s * ppcl->dN4_dy;
			de12_s = (pn1->dux_s * ppcl->dN1_dy + pn2->dux_s * ppcl->dN2_dy
					+ pn3->dux_s * ppcl->dN3_dy + pn4->dux_s * ppcl->dN4_dy
					+ pn1->duy_s * ppcl->dN1_dx + pn2->duy_s * ppcl->dN2_dx
					+ pn3->duy_s * ppcl->dN3_dx + pn4->duy_s * ppcl->dN4_dx) * 0.5;
			dw12 = (pn1->dux_s * ppcl->dN1_dy + pn2->dux_s * ppcl->dN2_dy
				  + pn3->dux_s * ppcl->dN3_dy + pn4->dux_s * ppcl->dN4_dy
				  - pn1->duy_s * ppcl->dN1_dx - pn2->duy_s * ppcl->dN2_dx
				  - pn3->duy_s * ppcl->dN3_dx - pn4->duy_s * ppcl->dN4_dx) * 0.5;

			// update strain (also assume that strain increment is Jaumann rate)
			//ppcl->de11 +=  dw12 * e12 * 2.0;
			//ppcl->de22 += -dw12 * e12 * 2.0;
			//ppcl->de12 +=  dw12 * (e22 - e11);
			ppcl->e11 += de11_s;
			ppcl->e22 += de22_s;
			ppcl->e12 += de12_s;

			// update stress
			E_tmp = ppcl->E / (1.0 + ppcl->niu) / (1.0 - 2.0 * ppcl->niu);
			ds11 = E_tmp * ((1.0 - ppcl->niu) * de11_s + ppcl->niu * de22_s);
			ds22 = E_tmp * (ppcl->niu * de11_s + (1.0 - ppcl->niu) * de22_s);
			ds12 = 2.0 * de12_s * ppcl->E / (2.0 * (1.0 + ppcl->niu));
			
			/* ------------------------------------------------------------------
			Rotate as Jaumann rate:
				tensor_rate = tensor_Jaumann_rate + tensor * dW_T + dW * tensor
			  ------------------------------------------------------------------- */
			/*			  ds11 +=  ppcl->dw12 * ppcl->s12 * 2.0;
						  ds22 += -ppcl->dw12 * ppcl->s12 * 2.0;
						  ds12 +=  ppcl->dw12 * (ppcl->s22 - ppcl->s11);	*/
			ppcl->s11 += ds11;
			ppcl->s22 += ds22;
			ppcl->s12 += ds12;

			// volumetric strain of solid phase
			de_vol_s = de11_s + de22_s;
			// volume
			ppcl->vol *= (1.0 + de_vol_s);
			// porosity
			ppcl->n = (de_vol_s + ppcl->n) / (1.0 + de_vol_s);
			// "volumetric strain" of fluid phase
			de_vol_f = -(1.0 - ppcl->n) / ppcl->n * de_vol_s
				   - (pn1->dux_f * ppcl->dN1_dx + pn2->dux_f * ppcl->dN2_dx + pn3->dux_f * ppcl->dN3_dx + pn4->dux_f * ppcl->dN4_dx)
				   - (pn1->duy_f * ppcl->dN1_dy + pn2->duy_f * ppcl->dN2_dy + pn3->duy_f * ppcl->dN3_dy + pn4->duy_f * ppcl->dN4_dy);

			// pore pressure
			ppcl->p += ppcl->Kf * de_vol_f;
			// fluid density
			ppcl->density_f += ppcl->density_f * de_vol_f;
		}
	}

	return 0;
}