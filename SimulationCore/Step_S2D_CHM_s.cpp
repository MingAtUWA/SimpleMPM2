#include "SimulationCore_pcp.h"

#include <cmath>
#include "Step_S2D_CHM_s.h"

namespace
{
	typedef Model_S2D_CHM_s::Particle Particle_mpm;
	typedef Model_S2D_CHM_s::ParticleCalVar ParticleCalVar_mpm;
	typedef Model_S2D_CHM_s::Element Element_mpm;
	typedef Model_S2D_CHM_s::Node Node_mpm;
};

Step_S2D_CHM_s::Step_S2D_CHM_s() :
	Step(&solve_substep_S2D_CHM_s),
	h_elem_raio(0.05), h_pcl_ratio(0.1),
	model(nullptr) {}

Step_S2D_CHM_s::~Step_S2D_CHM_s() {}

int Step_S2D_CHM_s::init_calculation(void)
{
	if (is_first_step)
	{
		for (size_t i = 0; i < model->pcl_num; ++i)
		{
			model->pcls[i].elem_num = 1;
		}
	}

	for (size_t i = 0; i < model->pcl_num; ++i)
	{
		Particle_mpm &pcl = model->pcls[i];
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux_s = 0.0;
		pcl.uy_s = 0.0;
		pcl.ux_f = 0.0;
		pcl.uy_f = 0.0;
	}
	return 0;
}

int Step_S2D_CHM_s::finalize_calculation(void) { return 0; }

int solve_substep_S2D_CHM_s(void *_self)
{
	Step_S2D_CHM_s &self = *(Step_S2D_CHM_s *)(_self);
	Model_S2D_CHM_s &model = *self.model;

	// init nodes
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Node_mpm &n = model.nodes[n_id];
		// solid phase
		n.m_s = 0.0;
		n.vx_s = 0.0;
		n.vy_s = 0.0;
		n.fx_kin_f = 0.0;
		n.fy_kin_f = 0.0;
		n.fx_ext_m = 0.0;
		n.fy_ext_m = 0.0;
		n.fx_int_m = 0.0;
		n.fy_int_m = 0.0;
		// fluid phase
		n.m_tf = 0.0;
		n.vx_f = 0.0;
		n.vy_f = 0.0;
		n.fx_ext_tf = 0.0;
		n.fy_ext_tf = 0.0;
		n.fx_int_tf = 0.0;
		n.fy_int_tf = 0.0;
		n.fx_drag_tf = 0.0;
		n.fy_drag_tf = 0.0;
	}

	// init particles
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = model.pcls[pcl_id];
		ParticleCalVar_mpm &pcl_var = pcl.var;
		if (pcl.elem_num)
		{
			if (!model.init_pcl_standard(pcl))
				continue;

			double m_s = (1.0 - pcl.n) * pcl.density_s * pcl_var.vol;
			double mvx_s = m_s * pcl.vx_s;
			double mvy_s = m_s * pcl.vy_s;
			double m_tf = pcl.density_f * pcl_var.vol;
			double mvx_tf = m_tf * pcl.vx_f;
			double mvy_tf = m_tf * pcl.vy_f;
			double n_prod_miu_div_k = pcl.n * pcl.miu / pcl.k;
			double n_miu_div_k_vrx_vol = n_prod_miu_div_k * (pcl.vx_f - pcl.vx_s) * pcl_var.vol;
			double n_miu_div_k_vry_vol = n_prod_miu_div_k * (pcl.vy_f - pcl.vy_s) * pcl_var.vol;

			// node 1
			Node_mpm &n1 = *pcl_var.pn1;
			// mixture phase
			n1.m_s += pcl_var.N1 * m_s;
			n1.vx_s += pcl_var.N1 * mvx_s;
			n1.vy_s += pcl_var.N1 * mvy_s;
			n1.fx_int_m += (pcl_var.dN1_dx * (pcl.s11 - pcl.p) + pcl_var.dN1_dy * pcl.s12) * pcl_var.vol;
			n1.fy_int_m += (pcl_var.dN1_dx * pcl.s12 + pcl_var.dN1_dy * (pcl.s22 - pcl.p)) * pcl_var.vol;
			// fluid phase
			n1.m_tf += pcl_var.N1 * m_tf;
			n1.vx_f += pcl_var.N1 * mvx_tf;
			n1.vy_f += pcl_var.N1 * mvy_tf;
			n1.fx_int_tf += (pcl_var.dN1_dx * -pcl.p) * pcl_var.vol;
			n1.fy_int_tf += (pcl_var.dN1_dy * -pcl.p) * pcl_var.vol;
			n1.fx_drag_tf += pcl_var.N1 * n_miu_div_k_vrx_vol;
			n1.fy_drag_tf += pcl_var.N1 * n_miu_div_k_vry_vol;
			
			// node 2
			Node_mpm &n2 = *pcl_var.pn2;
			// mixture phase
			n2.m_s += pcl_var.N2 * m_s;
			n2.vx_s += pcl_var.N2 * mvx_s;
			n2.vy_s += pcl_var.N2 * mvy_s;
			n2.fx_int_m += (pcl_var.dN2_dx * (pcl.s11 - pcl.p) + pcl_var.dN2_dy * pcl.s12) * pcl_var.vol;
			n2.fy_int_m += (pcl_var.dN2_dx * pcl.s12 + pcl_var.dN2_dy * (pcl.s22 - pcl.p)) * pcl_var.vol;
			// fluid phase
			n2.m_tf += pcl_var.N2 * m_tf;
			n2.vx_f += pcl_var.N2 * mvx_tf;
			n2.vy_f += pcl_var.N2 * mvy_tf;
			n2.fx_int_tf += (pcl_var.dN2_dx * -pcl.p) * pcl_var.vol;
			n2.fy_int_tf += (pcl_var.dN2_dy * -pcl.p) * pcl_var.vol;
			n2.fx_drag_tf += pcl_var.N2 * n_miu_div_k_vrx_vol;
			n2.fy_drag_tf += pcl_var.N2 * n_miu_div_k_vry_vol;

			// node 3
			Node_mpm &n3 = *pcl_var.pn3;
			// mixture phase
			n3.m_s += pcl_var.N3 * m_s;
			n3.vx_s += pcl_var.N3 * mvx_s;
			n3.vy_s += pcl_var.N3 * mvy_s;
			n3.fx_int_m += (pcl_var.dN3_dx * (pcl.s11 - pcl.p) + pcl_var.dN3_dy * pcl.s12) * pcl_var.vol;
			n3.fy_int_m += (pcl_var.dN3_dx * pcl.s12 + pcl_var.dN3_dy * (pcl.s22 - pcl.p)) * pcl_var.vol;
			// fluid phase
			n3.m_tf += pcl_var.N3 * m_tf;
			n3.vx_f += pcl_var.N3 * mvx_tf;
			n3.vy_f += pcl_var.N3 * mvy_tf;
			n3.fx_int_tf += (pcl_var.dN3_dx * -pcl.p) * pcl_var.vol;
			n3.fy_int_tf += (pcl_var.dN3_dy * -pcl.p) * pcl_var.vol;
			n3.fx_drag_tf += pcl_var.N3 * n_miu_div_k_vrx_vol;
			n3.fy_drag_tf += pcl_var.N3 * n_miu_div_k_vry_vol;

			// ------------------- node 4 -------------------
			Node_mpm &n4 = *pcl_var.pn4;
			// mixture phase
			n4.m_s += pcl_var.N4 * m_s;
			n4.vx_s += pcl_var.N4 * mvx_s;
			n4.vy_s += pcl_var.N4 * mvy_s;
			n4.fx_int_m += (pcl_var.dN4_dx * (pcl.s11 - pcl.p) + pcl_var.dN4_dy * pcl.s12) * pcl_var.vol;
			n4.fy_int_m += (pcl_var.dN4_dx * pcl.s12 + pcl_var.dN4_dy * (pcl.s22 - pcl.p)) * pcl_var.vol;
			// fluid phase
			n4.m_tf += pcl_var.N4 * m_tf;
			n4.vx_f += pcl_var.N4 * mvx_tf;
			n4.vy_f += pcl_var.N4 * mvy_tf;
			n4.fx_int_tf += (pcl_var.dN4_dx * -pcl.p) * pcl_var.vol;
			n4.fy_int_tf += (pcl_var.dN4_dy * -pcl.p) * pcl_var.vol;
			n4.fx_drag_tf += pcl_var.N4 * n_miu_div_k_vrx_vol;
			n4.fy_drag_tf += pcl_var.N4 * n_miu_div_k_vry_vol;
		}
	}

	// body force
	double bf_m, bf_tf;
	for (size_t bf_id = 0; bf_id < model.bfx_num; ++bf_id)
	{
		BodyForce &bf = model.bfxs[bf_id];
		Particle_mpm &pcl = model.pcls[bf.pcl_id];
		if (pcl.elem_num)
		{
			ParticleCalVar_mpm &pcl_var = pcl.var;
			// body force on particle
			bf_m = ((1.0 - pcl.n) * pcl.density_s + pcl.n * pcl.density_f) * pcl_var.vol * bf.bf;
			bf_tf = pcl.density_f * pcl_var.vol * bf.bf;
			// node 1
			Node_mpm &n1 = *pcl_var.pn1;
			n1.fx_ext_m += pcl_var.N1 * bf_m;
			n1.fx_ext_tf += pcl_var.N1 * bf_tf;
			// node 2
			Node_mpm &n2 = *pcl_var.pn2;
			n2.fx_ext_m += pcl_var.N2 * bf_m;
			n2.fx_ext_tf += pcl_var.N2 * bf_tf;
			// node 3
			Node_mpm &n3 = *pcl_var.pn3;
			n3.fx_ext_m += pcl_var.N3 * bf_m;
			n3.fx_ext_tf += pcl_var.N3 * bf_tf;
			// node 4
			Node_mpm &n4 = *pcl_var.pn4;
			n4.fx_ext_m += pcl_var.N4 * bf_m;
			n4.fx_ext_tf += pcl_var.N4 * bf_tf;
		}
	}
	for (size_t bf_id = 0; bf_id < model.bfy_num; ++bf_id)
	{
		BodyForce &bf = model.bfys[bf_id];
		Particle_mpm &pcl = model.pcls[bf.pcl_id];
		if (pcl.elem_num)
		{
			ParticleCalVar_mpm &pcl_var = pcl.var;
			// body force on particle
			bf_m = ((1.0 - pcl.n) * pcl.density_s + pcl.n * pcl.density_f) * pcl_var.vol * bf.bf;
			bf_tf = pcl.density_f * pcl_var.vol * bf.bf;
			// node 1
			Node_mpm &n1 = *pcl_var.pn1;
			n1.fy_ext_m += pcl_var.N1 * bf_m;
			n1.fy_ext_tf += pcl_var.N1 * bf_tf;
			// node 2
			Node_mpm &n2 = *pcl_var.pn2;
			n2.fy_ext_m += pcl_var.N2 * bf_m;
			n2.fy_ext_tf += pcl_var.N2 * bf_tf;
			// node 3
			Node_mpm &n3 = *pcl_var.pn3;
			n3.fy_ext_m += pcl_var.N3 * bf_m;
			n3.fy_ext_tf += pcl_var.N3 * bf_tf;
			// node 4
			Node_mpm &n4 = *pcl_var.pn4;
			n4.fy_ext_m += pcl_var.N4 * bf_m;
			n4.fy_ext_tf += pcl_var.N4 * bf_tf;
		}
	}

	// surface force
	for (size_t tf_id = 0; tf_id < model.tx_num; ++tf_id)
	{
		TractionBC_MPM &tf = model.txs[tf_id];
		Particle_mpm &pcl = model.pcls[tf.pcl_id];
		if (pcl.elem_num)
		{
			ParticleCalVar_mpm &pcl_var = pcl.var;
			// node 1
			Node_mpm &n1 = *pcl_var.pn1;
			n1.fx_ext_m += pcl_var.N1 * tf.t;
			// node 2
			Node_mpm &n2 = *pcl_var.pn2;
			n2.fx_ext_m += pcl_var.N2 * tf.t;
			// node 3
			Node_mpm &n3 = *pcl_var.pn3;
			n3.fx_ext_m += pcl_var.N3 * tf.t;
			// node 4
			Node_mpm &n4 = *pcl_var.pn4;
			n4.fx_ext_m += pcl_var.N4 * tf.t;
		}
	}
	for (size_t tf_id = 0; tf_id < model.ty_num; ++tf_id)
	{
		TractionBC_MPM &tf = model.tys[tf_id];
		Particle_mpm &pcl = model.pcls[tf.pcl_id];
		if (pcl.elem_num)
		{
			ParticleCalVar_mpm &pcl_var = pcl.var;
			// node 1
			Node_mpm &n1 = *pcl_var.pn1;
			n1.fy_ext_m += pcl_var.N1 * tf.t;
			// node 2
			Node_mpm &n2 = *pcl_var.pn2;
			n2.fy_ext_m += pcl_var.N2 * tf.t;
			// node 3
			Node_mpm &n3 = *pcl_var.pn3;
			n3.fy_ext_m += pcl_var.N3 * tf.t;
			// node 4
			Node_mpm &n4 = *pcl_var.pn4;
			n4.fy_ext_m += pcl_var.N4 * tf.t;
		}
	}
	// pore pressure force...

	// update nodal acceleration of fluid pahse
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Node_mpm &n = model.nodes[n_id];
		if (n.m_tf != 0.0)
		{
			n.ax_f = (n.fx_ext_tf - n.fx_int_tf - n.fx_drag_tf) / n.m_tf;
			n.ay_f = (n.fy_ext_tf - n.fy_int_tf - n.fy_drag_tf) / n.m_tf;
		}
	}
	for (size_t a_id = 0; a_id < model.asx_num; ++a_id)
	{
		Node_mpm &n = model.nodes[model.afxs[a_id].node_id];
		n.ax_f = model.afxs[a_id].a;
	}
	for (size_t a_id = 0; a_id < model.asy_num; ++a_id)
	{
		Node_mpm &n = model.nodes[model.afys[a_id].node_id];
		n.ay_f = model.afys[a_id].a;
	}

	// update nodal momentum of fluid phase
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Node_mpm &n = model.nodes[n_id];
		if (n.m_tf != 0.0)
		{
			n.vx_f /= n.m_tf;
			n.vx_f += n.ax_f * self.dtime;
			n.vy_f /= n.m_tf;
			n.vy_f += n.ay_f * self.dtime;
		}
	}
	// apply velocity boundary conditions of fluid phase
	for (size_t n_id = 0; n_id < model.vfx_num; ++n_id)
	{
		Node_mpm &n =  model.nodes[model.vfxs[n_id].node_id];
		n.vx_f = model.vfxs[n_id].v;
		n.ax_f = 0.0;
	}
	for (size_t n_id = 0; n_id < model.vfy_num; ++n_id)
	{
		Node_mpm &n = model.nodes[model.vfys[n_id].node_id];
		n.vy_f = model.vfys[n_id].v;
		n.ay_f = 0.0;
	}

	// calculate the inertial term of fluid in mixture formulation
	double pcl_ax_f, pcl_ay_f;
	double pcl_m_f, pcl_max_f, pcl_may_f;
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = model.pcls[pcl_id];
		if (pcl.elem_num)
		{
			ParticleCalVar_mpm &pcl_var = pcl.var;
			Node_mpm &n1 = *pcl_var.pn1;
			Node_mpm &n2 = *pcl_var.pn2;
			Node_mpm &n3 = *pcl_var.pn3;
			Node_mpm &n4 = *pcl_var.pn4;
			// particle acceleration
			pcl_ax_f = pcl_var.N1 * n1.ax_f + pcl_var.N2 * n2.ax_f
					 + pcl_var.N3 * n3.ax_f + pcl_var.N4 * n4.ax_f;
			pcl_ay_f = pcl_var.N1 * n1.ay_f + pcl_var.N2 * n2.ay_f
					 + pcl_var.N3 * n3.ay_f + pcl_var.N4 * n4.ay_f;
			pcl_m_f = pcl.n * pcl.density_f * pcl_var.vol;
			pcl_max_f = pcl_m_f * pcl_ax_f;
			pcl_may_f = pcl_m_f * pcl_ay_f;
			// node 1
			n1.fx_kin_f += pcl_var.N1 * pcl_max_f;
			n1.fy_kin_f += pcl_var.N1 * pcl_may_f;
			// node 2
			n2.fx_kin_f += pcl_var.N2 * pcl_max_f;
			n2.fy_kin_f += pcl_var.N2 * pcl_may_f;
			// node 3
			n3.fx_kin_f += pcl_var.N3 * pcl_max_f;
			n3.fy_kin_f += pcl_var.N3 * pcl_may_f;
			// node 4
			n4.fx_kin_f += pcl_var.N4 * pcl_max_f;
			n4.fy_kin_f += pcl_var.N4 * pcl_may_f;
		}
	}

	// update nodal velocity of solid phase
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Node_mpm &n = model.nodes[n_id];
		if (n.m_s != 0.0)
		{
			n.ax_s = (n.fx_ext_m - n.fx_int_m - n.fx_kin_f) / n.m_s;
			n.ay_s = (n.fy_ext_m - n.fy_int_m - n.fy_kin_f) / n.m_s;
		}
	}
	// apply acceleration boundary conditions
	for (size_t a_id = 0; a_id < model.asx_num; ++a_id)
	{
		Node_mpm &n = model.nodes[model.asxs[a_id].node_id];
		n.ax_s = model.asxs[a_id].a;
	}
	for (size_t a_id = 0; a_id < model.asy_num; ++a_id)
	{
		Node_mpm &n = model.nodes[model.asys[a_id].node_id];
		n.ay_s = model.asys[a_id].a;
	}

	// update nodal momentum of fluid pahse
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Node_mpm &n = model.nodes[n_id];
		if (n.m_s != 0.0)
		{
			n.vx_s /= n.m_s;
			n.vx_s += n.ax_s * self.dtime;
			n.vy_s /= n.m_s;
			n.vy_s += n.ay_s * self.dtime;
		}
	}
	// apply velocity boundary conditions of solid phase
	for (size_t v_id = 0; v_id < model.vsx_num; ++v_id)
	{
		Node_mpm &n = model.nodes[model.vsxs[v_id].node_id];
		n.vx_s = model.vsxs[v_id].v;
		n.ax_s = 0.0;
	}
	for (size_t v_id = 0; v_id < model.vsy_num; ++v_id)
	{
		Node_mpm &n = model.nodes[model.vsys[v_id].node_id];
		n.vy_s = model.vsys[v_id].v;
		n.ay_s = 0.0;
	}

	// update displacement increment of both phases
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Node_mpm &n = model.nodes[n_id];
		if (n.m_s != 0.0)
		{
			// solid phase
			n.dux_s = n.vx_s * self.dtime;
			n.duy_s = n.vy_s * self.dtime;
			// fluid phase
			n.dux_f = n.vx_f * self.dtime;
			n.duy_f = n.vy_f * self.dtime;
		}
	}

	// map variables back to and update variables particles
	double de11, de22, de12, dw12;
	double ds11, ds22, ds12;
	double de_vol_s, de_vol_f;
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = model.pcls[pcl_id];
		if (pcl.elem_num)
		{
			ParticleCalVar_mpm &pcl_var = pcl.var;
			Node_mpm &n1 = *pcl_var.pn1;
			Node_mpm &n2 = *pcl_var.pn2;
			Node_mpm &n3 = *pcl_var.pn3;
			Node_mpm &n4 = *pcl_var.pn4;
			
			// velocity
			pcl.vx_s += (n1.ax_s * pcl_var.N1 + n2.ax_s * pcl_var.N2
					   + n3.ax_s * pcl_var.N3 + n4.ax_s * pcl_var.N4) * self.dtime;
			pcl.vy_s += (n1.ay_s * pcl_var.N1 + n2.ay_s * pcl_var.N2
					   + n3.ay_s * pcl_var.N3 + n4.ay_s * pcl_var.N4) * self.dtime;
			pcl.vx_f += (n1.ax_f * pcl_var.N1 + n2.ax_f * pcl_var.N2
					   + n3.ax_f * pcl_var.N3 + n4.ax_f * pcl_var.N4) * self.dtime;
			pcl.vy_f += (n1.ay_f * pcl_var.N1 + n2.ay_f * pcl_var.N2
					   + n3.ay_f * pcl_var.N3 + n4.ay_f * pcl_var.N4) * self.dtime;

			// displacement
			pcl.ux_s += n1.dux_s * pcl_var.N1 + n2.dux_s * pcl_var.N2
					  + n3.dux_s * pcl_var.N3 + n4.dux_s * pcl_var.N4;
			pcl.uy_s += n1.duy_s * pcl_var.N1 + n2.duy_s * pcl_var.N2
					  + n3.duy_s * pcl_var.N3 + n4.duy_s * pcl_var.N4;
			pcl.ux_f += n1.dux_f * pcl_var.N1 + n2.dux_f * pcl_var.N2
					  + n3.dux_f * pcl_var.N3 + n4.dux_f * pcl_var.N4;
			pcl.uy_f += n1.duy_f * pcl_var.N1 + n2.duy_f * pcl_var.N2
					  + n3.duy_f * pcl_var.N3 + n4.duy_f * pcl_var.N4;

			// update position
			pcl.x = pcl.x_ori + pcl.ux_s;
			pcl.y = pcl.y_ori + pcl.uy_s;

			// strain increment
			de11 = n1.dux_s * pcl_var.dN1_dx + n2.dux_s * pcl_var.dN2_dx
				 + n3.dux_s * pcl_var.dN3_dx + n4.dux_s * pcl_var.dN4_dx;
			de22 = n1.duy_s * pcl_var.dN1_dy + n2.duy_s * pcl_var.dN2_dy
				 + n3.duy_s * pcl_var.dN3_dy + n4.duy_s * pcl_var.dN4_dy;
			de12 = (n1.dux_s * pcl_var.dN1_dy + n2.dux_s * pcl_var.dN2_dy
				  + n3.dux_s * pcl_var.dN3_dy + n4.dux_s * pcl_var.dN4_dy
				  + n1.duy_s * pcl_var.dN1_dx + n2.duy_s * pcl_var.dN2_dx
				  + n3.duy_s * pcl_var.dN3_dx + n4.duy_s * pcl_var.dN4_dx) * 0.5;
			dw12 = (n1.dux_s * pcl_var.dN1_dy + n2.dux_s * pcl_var.dN2_dy
				  + n3.dux_s * pcl_var.dN3_dy + n4.dux_s * pcl_var.dN4_dy
				  - n1.duy_s * pcl_var.dN1_dx - n2.duy_s * pcl_var.dN2_dx
				  - n3.duy_s * pcl_var.dN3_dx - n4.duy_s * pcl_var.dN4_dx) * 0.5;
			// update strain (also assume that strain increment is Jaumann rate)
			//de11 +=  dw12 * e12 * 2.0;
			//de22 += -dw12 * e12 * 2.0;
			//de12 +=  dw12 * (e22 - e11);
			pcl.e11 += de11;
			pcl.e22 += de22;
			pcl.e12 += de12;

			// update stress
			double E_tmp = pcl.E / (1.0 + pcl.niu) / (1.0 - 2.0 * pcl.niu);
			ds11 = E_tmp * ((1.0 - pcl.niu) * de11 + pcl.niu * de22);
			ds22 = E_tmp * (pcl.niu * de11 + (1.0 - pcl.niu) * de22);
			ds12 = 2.0 * pcl.E / (2.0 * (1.0 + pcl.niu)) * de12;
			/* ------------------------------------------------------------------
			Rotate as Jaumann rate:
			tensor_rate = tensor_Jaumann_rate + tensor * dW_T + dW * tensor
			------------------------------------------------------------------- */
			//ds11 +=  dw12 * pcl_var.s12 * 2.0;
			//ds22 += -dw12 * pcl_var.s12 * 2.0;
			//ds12 +=  dw12 * (pcl_var.s22 - pcl_var.s11);
			pcl.s11 += ds11;
			pcl.s22 += ds22;
			pcl.s12 += ds12;

			// volumetric strain of solid phase
			de_vol_s = de11 + de22;
			// "volumetric strain" of fluid phase
			de_vol_f = -(1.0 - pcl.n) / pcl.n * de_vol_s
					  - (n1.dux_f * pcl_var.dN1_dx + n2.dux_f * pcl_var.dN2_dx
					   + n3.dux_f * pcl_var.dN3_dx + n4.dux_f * pcl_var.dN4_dx)
					  - (n1.duy_f * pcl_var.dN1_dy + n2.duy_f * pcl_var.dN2_dy
					   + n3.duy_f * pcl_var.dN3_dy + n4.duy_f * pcl_var.dN4_dy);
			// porosity
			pcl.n = (de_vol_s + pcl.n) / (1.0 + de_vol_s);

			// fluid density
			pcl.density_f += pcl.density_f * de_vol_f;
			// pore pressure (can use EOS and get p from density_f instead)
			pcl.p += pcl.Kf * de_vol_f;
		}
	}
	
	return 0;
}
