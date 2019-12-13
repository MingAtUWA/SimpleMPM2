#include "SimulationCore_pcp.h"

#include <cmath>
#include "Step_T2D_CHM_s_SE.h"

Step_T2D_CHM_s_SE::Step_T2D_CHM_s_SE() :
	Step(&solve_substep_T2D_CHM_s_SE), model(nullptr),
	damping_ratio(0.0), bv_ratio(0.0) {}

Step_T2D_CHM_s_SE::~Step_T2D_CHM_s_SE() {}

int Step_T2D_CHM_s_SE::init_calculation(void)
{
	Model_T2D_CHM_s &md = *model;

	if (is_first_step) {}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		pcl.pe = md.elems;
		pcl.vol_s = pcl.m_s / pcl.density_s;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux_s = 0.0;
		pcl.uy_s = 0.0;
		pcl.ux_f = 0.0;
		pcl.uy_f = 0.0;
	}

	return 0;
}

int Step_T2D_CHM_s_SE::finalize_calculation(void) { return 0; }

namespace
{
	typedef Model_T2D_CHM_s::Particle Particle_mpm;
	typedef Model_T2D_CHM_s::Element Element_mpm;
	typedef Model_T2D_CHM_s::Node Node_mpm;
}

int solve_substep_T2D_CHM_s_SE(void *_self)
{
	Step_T2D_CHM_s_SE &self = *(Step_T2D_CHM_s_SE *)(_self);
	Model_T2D_CHM_s &md = *self.model;

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_mpm &n = md.nodes[n_id];
		// material point
		n.has_mp = false;
		// solid phase
		n.m_s = 0.0;
		n.vx_s = 0.0;
		n.vy_s = 0.0;
		n.fx_ext_s = 0.0;
		n.fy_ext_s = 0.0;
		n.fx_int_s = 0.0;
		n.fy_int_s = 0.0;
		// fluid phase
		n.m_f = 0.0;
		n.vx_f = 0.0;
		n.vy_f = 0.0;
		n.fx_ext_f = 0.0;
		n.fy_ext_f = 0.0;
		n.fx_int_f = 0.0;
		n.fy_int_f = 0.0;
		// solid - fluid interaction
		n.fx_drag = 0.0;
		n.fy_drag = 0.0;
		// rigid body
		n.has_rb = false;
		n.vol_rb = 0.0;
		n.vx_rb = 0.0;
		n.vy_rb = 0.0;
		// strain enhancement approach
		n.vol = 0.0;
		n.de_vol_s = 0.0;
		n.de_vol_f = 0.0;
	}

	// init elements
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element_mpm &e = md.elems[e_id];
		e.vol = 0.0;
		e.n = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
		e.p = 0.0;
		e.pcls = nullptr;
	}

	// init particles
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			if (!(pcl.pe = md.find_in_which_element(pcl)))
				continue;
			pcl.pe->add_pcl(pcl);

			double mvx_s = pcl.m_s * pcl.vx_s;
			double mvy_s = pcl.m_s * pcl.vy_s;
			pcl.vol = pcl.vol_s / (1.0 - pcl.n);
			pcl.m_f = pcl.n * pcl.density_f * pcl.vol;
			double mvx_f = pcl.m_f * pcl.vx_f;
			double mvy_f = pcl.m_f * pcl.vy_f;
			double n2_miu_div_k = pcl.n * pcl.n * md.miu / md.k;
			double n2_miu_div_k_vrx_vol = n2_miu_div_k * (pcl.vx_f - pcl.vx_s) * pcl.vol;
			double n2_miu_div_k_vry_vol = n2_miu_div_k * (pcl.vy_f - pcl.vy_s) * pcl.vol;
			
			Element_mpm &e = *pcl.pe;
			e.vol += pcl.vol;
			e.n += pcl.vol_s;
			e.s11 += pcl.vol * pcl.s11;
			e.s22 += pcl.vol * pcl.s22;
			e.s12 += pcl.vol * pcl.s12;
			e.p += pcl.vol * pcl.p;

			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.has_mp = true;
			// solid phase
			n1.m_s += pcl.N1 * pcl.m_s;
			n1.vx_s += pcl.N1 * mvx_s;
			n1.vy_s += pcl.N1 * mvy_s;
			// fluid phase
			n1.m_f += pcl.N1 * pcl.m_f;
			n1.vx_f += pcl.N1 * mvx_f;
			n1.vy_f += pcl.N1 * mvy_f;
			// solid - fluid interaction
			n1.fx_drag += pcl.N1 * n2_miu_div_k_vrx_vol;
			n1.fy_drag += pcl.N1 * n2_miu_div_k_vry_vol;

			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.has_mp = true;
			// mixture phase
			n2.m_s += pcl.N2 * pcl.m_s;
			n2.vx_s += pcl.N2 * mvx_s;
			n2.vy_s += pcl.N2 * mvy_s;
			// fluid phase
			n2.m_f += pcl.N2 * pcl.m_f;
			n2.vx_f += pcl.N2 * mvx_f;
			n2.vy_f += pcl.N2 * mvy_f;
			// solid - fluid interaction
			n2.fx_drag += pcl.N2 * n2_miu_div_k_vrx_vol;
			n2.fy_drag += pcl.N2 * n2_miu_div_k_vry_vol;

			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.has_mp = true;
			// mixture phase
			n3.m_s += pcl.N3 * pcl.m_s;
			n3.vx_s += pcl.N3 * mvx_s;
			n3.vy_s += pcl.N3 * mvy_s;
			// fluid phase
			n3.m_f += pcl.N3 * pcl.m_f;
			n3.vx_f += pcl.N3 * mvx_f;
			n3.vy_f += pcl.N3 * mvy_f;
			// solid - fluid interaction
			n3.fx_drag += pcl.N3 * n2_miu_div_k_vrx_vol;
			n3.fy_drag += pcl.N3 * n2_miu_div_k_vry_vol;
		}
	}

	// update nodal momentum
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_mpm &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.vx_s /= n.m_s;
			n.vy_s /= n.m_s;
			n.vx_f /= n.m_f;
			n.vy_f /= n.m_f;
		}
	}

	// apply velocity bc
	for (size_t v_id = 0; v_id < md.vsx_num; ++v_id)
	{
		Node_mpm &n = md.nodes[md.vsxs[v_id].node_id];
		n.vx_s = md.vsxs[v_id].v;
	}
	for (size_t v_id = 0; v_id < md.vsy_num; ++v_id)
	{
		Node_mpm &n = md.nodes[md.vsys[v_id].node_id];
		n.vy_s = md.vsys[v_id].v;
	}
	for (size_t v_id = 0; v_id < md.vfx_num; ++v_id)
	{
		Node_mpm &n = md.nodes[md.vfxs[v_id].node_id];
		n.vx_f = md.vfxs[v_id].v;
	}
	for (size_t v_id = 0; v_id < md.vfy_num; ++v_id)
	{
		Node_mpm &n = md.nodes[md.vfys[v_id].node_id];
		n.vy_f = md.vfys[v_id].v;
	}

	double de_vol_f_rate;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element_mpm &e = md.elems[e_id];
		if (e.pcls && e.vol != 0.0)
		{
			e.n = 1.0 - e.n / e.vol; // 1.0 - Vs / V
			e.s11 /= e.vol;
			e.s22 /= e.vol;
			e.s12 /= e.vol;
			e.p /= e.vol;
			if (e.vol > e.area)
				e.vol = e.area;

			Node_mpm &n1 = md.nodes[e.n1];
			Node_mpm &n2 = md.nodes[e.n2];
			Node_mpm &n3 = md.nodes[e.n3];

			//// volume strain rate
			//de_vol_f_rate = -(1.0 - e.n) / e.n *
			//				 (n1.vx_s * e.dN1_dx + n2.vx_s * e.dN2_dx + n3.vx_s * e.dN3_dx
			//				+ n1.vy_s * e.dN1_dy + n2.vy_s * e.dN2_dy + n3.vy_s * e.dN3_dy)
			//				-(n1.vx_f * e.dN1_dx + n2.vx_f * e.dN2_dx + n3.vx_f * e.dN3_dx
			//				+ n1.vy_f * e.dN1_dy + n2.vy_f * e.dN2_dy + n3.vy_f * e.dN3_dy);
			//// add bulk viscosity
			//e.p -= self.bv_ratio * de_vol_f_rate;
			
			// node 1
			n1.fx_int_s += (e.dN1_dx * (e.s11 - (1.0 - e.n) * e.p) + e.dN1_dy * e.s12) * e.vol;
			n1.fy_int_s += (e.dN1_dx * e.s12 + e.dN1_dy * (e.s22 - (1.0 - e.n) * e.p)) * e.vol;
			n1.fx_int_f += (e.dN1_dx * e.n * -e.p) * e.vol;
			n1.fy_int_f += (e.dN1_dy * e.n * -e.p) * e.vol;
			// node 2
			n2.fx_int_s += (e.dN2_dx * (e.s11 - (1.0 - e.n) * e.p) + e.dN2_dy * e.s12) * e.vol;
			n2.fy_int_s += (e.dN2_dx * e.s12 + e.dN2_dy * (e.s22 - (1.0 - e.n) * e.p)) * e.vol;
			n2.fx_int_f += (e.dN2_dx * e.n * -e.p) * e.vol;
			n2.fy_int_f += (e.dN2_dy * e.n * -e.p) * e.vol;
			// node 3
			n3.fx_int_s += (e.dN3_dx * (e.s11 - (1.0 - e.n) * e.p) + e.dN3_dy * e.s12) * e.vol;
			n3.fy_int_s += (e.dN3_dx * e.s12 + e.dN3_dy * (e.s22 - (1.0 - e.n) * e.p)) * e.vol;
			n3.fx_int_f += (e.dN3_dx * e.n * -e.p) * e.vol;
			n3.fy_int_f += (e.dN3_dy * e.n * -e.p) * e.vol;
		}
	}

	// body force
	double bf_s, bf_f;
	for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
	{
		BodyForce &bf = md.bfxs[bf_id];
		Particle_mpm &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			// body force on particle
			bf_s = pcl.m_s * bf.bf;
			bf_f = pcl.m_f * bf.bf;
			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.fx_ext_s += pcl.N1 * bf_s;
			n1.fx_ext_f += pcl.N1 * bf_f;
			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.fx_ext_s += pcl.N2 * bf_s;
			n2.fx_ext_f += pcl.N2 * bf_f;
			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.fx_ext_s += pcl.N3 * bf_s;
			n3.fx_ext_f += pcl.N3 * bf_f;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
	{
		BodyForce &bf = md.bfys[bf_id];
		Particle_mpm &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			// body force on particle
			bf_s = pcl.m_s * bf.bf;
			bf_f = pcl.m_f * bf.bf;
			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.fy_ext_s += pcl.N1 * bf_s;
			n1.fy_ext_f += pcl.N1 * bf_f;
			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.fy_ext_s += pcl.N2 * bf_s;
			n2.fy_ext_f += pcl.N2 * bf_f;
			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.fy_ext_s += pcl.N3 * bf_s;
			n3.fy_ext_f += pcl.N3 * bf_f;
		}
	}

	// surface force
	for (size_t tf_id = 0; tf_id < md.tx_num; ++tf_id)
	{
		TractionBC_MPM &tf = md.txs[tf_id];
		Particle_mpm &pcl = md.pcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.fx_ext_s += pcl.N1 * tf.t;
			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.fx_ext_s += pcl.N2 * tf.t;
			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.fx_ext_s += pcl.N3 * tf.t;
		}
	}
	for (size_t tf_id = 0; tf_id < md.ty_num; ++tf_id)
	{
		TractionBC_MPM &tf = md.tys[tf_id];
		Particle_mpm &pcl = md.pcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.fy_ext_s += pcl.N1 * tf.t;
			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.fy_ext_s += pcl.N2 * tf.t;
			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.fy_ext_s += pcl.N3 * tf.t;
		}
	}
	// pore pressure force...

	// update nodal acceleration of fluid pahse
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_mpm &n = md.nodes[n_id];
		double nf, v_sign;
		if (n.has_mp) // or n.m_f != 0.0
		{
			// fx_s
			if (n.vx_s > 0.0)
				v_sign = 1.0;
			else if (n.vx_s < 0.0)
				v_sign = -1.0;
			else
				v_sign = 0.0;
			nf = n.fx_ext_s - n.fx_int_s;
			n.ax_s = (nf + n.fx_drag - self.damping_ratio * abs(nf) * v_sign) / n.m_s;
			// fy_s
			if (n.vy_s > 0.0)
				v_sign = 1.0;
			else if (n.vy_s < 0.0)
				v_sign = -1.0;
			else
				v_sign = 0.0;
			nf = n.fy_ext_s - n.fy_int_s;
			n.ay_s = (nf + n.fy_drag - self.damping_ratio * abs(nf) * v_sign) / n.m_s;
			// fx_f
			if (n.vx_f > 0.0)
				v_sign = 1.0;
			else if (n.vx_f < 0.0)
				v_sign = -1.0;
			else
				v_sign = 0.0;
			nf = n.fx_ext_f - n.fx_int_f;
			n.ax_f = (nf - n.fx_drag - self.damping_ratio * abs(nf) * v_sign) / n.m_f;
			// fy_f
			if (n.vy_f > 0.0)
				v_sign = 1.0;
			else if (n.vy_f < 0.0)
				v_sign = -1.0;
			else
				v_sign = 0.0;
			nf = n.fy_ext_f - n.fy_int_f;
			n.ay_f = (nf - n.fy_drag - self.damping_ratio * abs(nf) * v_sign) / n.m_f;
		}
	}

	// apply acceleration bc
	for (size_t a_id = 0; a_id < md.asx_num; ++a_id)
	{
		Node_mpm &n = md.nodes[md.asxs[a_id].node_id];
		n.ax_s = md.asxs[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.asy_num; ++a_id)
	{
		Node_mpm &n = md.nodes[md.asys[a_id].node_id];
		n.ay_s = md.asys[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.afx_num; ++a_id)
	{
		Node_mpm &n = md.nodes[md.afxs[a_id].node_id];
		n.ax_f = md.afxs[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.afy_num; ++a_id)
	{
		Node_mpm &n = md.nodes[md.afys[a_id].node_id];
		n.ay_f = md.afys[a_id].a;
	}

	// update nodal momentum
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_mpm &n = md.nodes[n_id];
		if (n.has_mp && n.m_s != 0.0)
		{
			n.vx_s += n.ax_s * self.dtime;
			n.vy_s += n.ay_s * self.dtime;
			n.vx_f += n.ax_f * self.dtime;
			n.vy_f += n.ay_f * self.dtime;
		}
	}

	// contact detection and velocity modification
	if (md.get_rigid_circle().is_init())
		md.apply_rigid_body_to_bg_mesh(self.dtime);

	// apply velocity bc
	for (size_t v_id = 0; v_id < md.vsx_num; ++v_id)
	{
		Node_mpm &n = md.nodes[md.vsxs[v_id].node_id];
		n.vx_s = md.vsxs[v_id].v;
		n.ax_s = 0.0;
	}
	for (size_t v_id = 0; v_id < md.vsy_num; ++v_id)
	{
		Node_mpm &n = md.nodes[md.vsys[v_id].node_id];
		n.vy_s = md.vsys[v_id].v;
		n.ay_s = 0.0;
	}
	for (size_t v_id = 0; v_id < md.vfx_num; ++v_id)
	{
		Node_mpm &n =  md.nodes[md.vfxs[v_id].node_id];
		n.vx_f = md.vfxs[v_id].v;
		n.ax_f = 0.0;
	}
	for (size_t v_id = 0; v_id < md.vfy_num; ++v_id)
	{
		Node_mpm &n = md.nodes[md.vfys[v_id].node_id];
		n.vy_f = md.vfys[v_id].v;
		n.ay_f = 0.0;
	}

	// update displacement increment of both phases
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_mpm &n = md.nodes[n_id];
		if (n.has_mp && n.m_s != 0.0)
		{
			// solid phase
			n.dux_s = n.vx_s * self.dtime;
			n.duy_s = n.vy_s * self.dtime;
			// fluid phase
			n.dux_f = n.vx_f * self.dtime;
			n.duy_f = n.vy_f * self.dtime;
		}
	}

	// map variables back to particles and update their variables
	double de11, de22, de_mean;
	double ds11, ds22, ds12;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element_mpm &e = md.elems[e_id];
		if (e.pcls)
		{
			Node_mpm &n1 = md.nodes[e.n1];
			Node_mpm &n2 = md.nodes[e.n2];
			Node_mpm &n3 = md.nodes[e.n3];

			// strain increment
			de11 = n1.dux_s * e.dN1_dx + n2.dux_s * e.dN2_dx + n3.dux_s * e.dN3_dx;
			de22 = n1.duy_s * e.dN1_dy + n2.duy_s * e.dN2_dy + n3.duy_s * e.dN3_dy;
			e.de12 = (n1.dux_s * e.dN1_dy + n2.dux_s * e.dN2_dy + n3.dux_s * e.dN3_dy
			 		+ n1.duy_s * e.dN1_dx + n2.duy_s * e.dN2_dx + n3.duy_s * e.dN3_dx) * 0.5;
			// volumetric strain of solid phase
			e.de_vol_s = de11 + de22;
			de_mean = e.de_vol_s / 3.0;
			e.dde11 = de11 - de_mean;
			e.dde22 = de22 - de_mean;
			// "volumetric strain" of fluid phase, take compression as positive
			e.de_vol_f = (1.0 - e.n) / e.n * -e.de_vol_s
					   - (n1.dux_f * e.dN1_dx + n2.dux_f * e.dN2_dx + n3.dux_f * e.dN3_dx
						+ n1.duy_f * e.dN1_dy + n2.duy_f * e.dN2_dy + n3.duy_f * e.dN3_dy);
		}
	}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			double vol_de_vol_s = pcl.vol * e.de_vol_s;
			double vol_de_vol_f = pcl.vol * e.de_vol_f;
			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.vol += pcl.N1 * pcl.vol;
			n1.de_vol_s += pcl.N1 * vol_de_vol_s;
			n1.de_vol_f += pcl.N1 * vol_de_vol_f;
			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.vol += pcl.N2 * pcl.vol;
			n2.de_vol_s += pcl.N2 * vol_de_vol_s;
			n2.de_vol_f += pcl.N2 * vol_de_vol_f;
			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.vol += pcl.N3 * pcl.vol;
			n3.de_vol_s += pcl.N3 * vol_de_vol_s;
			n3.de_vol_f += pcl.N3 * vol_de_vol_f;
		}
	}

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_mpm &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.de_vol_s /= n.vol;
			n.de_vol_f /= n.vol;
		}
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element_mpm &e = md.elems[e_id];
		if (e.pcls)
		{
			Node_mpm &n1 = md.nodes[e.n1];
			Node_mpm &n2 = md.nodes[e.n2];
			Node_mpm &n3 = md.nodes[e.n3];
			//e.de_vol_s = (n1.de_vol_s + n2.de_vol_s + n3.de_vol_s) / 3.0;
			e.de_vol_f = (n1.de_vol_f + n2.de_vol_f + n3.de_vol_f) / 3.0;
			e.p += md.Kf * e.de_vol_f;
		}
	}

	double de12;
	double de_vol_s, de_vol_f;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			Node_mpm &n1 = md.nodes[e.n1];
			Node_mpm &n2 = md.nodes[e.n2];
			Node_mpm &n3 = md.nodes[e.n3];

			// velocity
			//pcl.vx_s = n1.vx_s * pcl.N1 + n2.vx_s * pcl.N2 + n3.vx_s * pcl.N3;
			//pcl.vy_s = n1.vy_s * pcl.N1 + n2.vy_s * pcl.N2 + n3.vy_s * pcl.N3;
			//pcl.vx_f = n1.vx_f * pcl.N1 + n2.vx_f * pcl.N2 + n3.vx_f * pcl.N3;
			//pcl.vy_f = n1.vy_f * pcl.N1 + n2.vy_f * pcl.N2 + n3.vy_f * pcl.N3;
			pcl.vx_s += (n1.ax_s * pcl.N1 + n2.ax_s * pcl.N2 + n3.ax_s * pcl.N3) * self.dtime;
			pcl.vy_s += (n1.ay_s * pcl.N1 + n2.ay_s * pcl.N2 + n3.ay_s * pcl.N3) * self.dtime;
			pcl.vx_f += (n1.ax_f * pcl.N1 + n2.ax_f * pcl.N2 + n3.ax_f * pcl.N3) * self.dtime;
			pcl.vy_f += (n1.ay_f * pcl.N1 + n2.ay_f * pcl.N2 + n3.ay_f * pcl.N3) * self.dtime;

			// displacement
			pcl.ux_s += n1.dux_s * pcl.N1 + n2.dux_s * pcl.N2 + n3.dux_s * pcl.N3;
			pcl.uy_s += n1.duy_s * pcl.N1 + n2.duy_s * pcl.N2 + n3.duy_s * pcl.N3;
			pcl.ux_f += n1.dux_f * pcl.N1 + n2.dux_f * pcl.N2 + n3.dux_f * pcl.N3;
			pcl.uy_f += n1.duy_f * pcl.N1 + n2.duy_f * pcl.N2 + n3.duy_f * pcl.N3;

			// update position
			pcl.x = pcl.x_ori + pcl.ux_s;
			pcl.y = pcl.y_ori + pcl.uy_s;

			// strain enhancement appraoch
			de_vol_s = n1.de_vol_s * pcl.N1 + n2.de_vol_s * pcl.N2 + n3.de_vol_s * pcl.N3;
			// porosity
			pcl.n = (de_vol_s + pcl.n) / (1.0 + de_vol_s);

			// strain
			de_vol_s /= 3.0;
			de11 = e.dde11 + de_vol_s;
			de22 = e.dde22 + de_vol_s;
			de12 = e.de12;
			// strain
			pcl.e11 += de11;
			pcl.e22 += de22;
			pcl.e12 += de12;

			// update stress
			// need remapping back to yield surface
			double E_tmp = md.E / ((1.0 + md.niu) * (1.0 - 2.0 * md.niu));
			ds11 = E_tmp * ((1.0 - md.niu) * de11 + md.niu * de22);
			ds22 = E_tmp * (md.niu * de11 + (1.0 - md.niu) * de22);
			ds12 = md.E / (2.0 * (1.0 + md.niu)) * 2.0 * de12;
			pcl.s11 += ds11;
			pcl.s22 += ds22;
			pcl.s12 += ds12;

			// pore pressure
			//de_vol_f = n1.de_vol_f * pcl.N1 + n2.de_vol_f * pcl.N2 + n3.de_vol_f * pcl.N3;
			//pcl.p += md.Kf * de_vol_f;
			pcl.p = e.p;
			// fluid density
			pcl.density_f /= 1.0 - e.de_vol_f;
			//pcl.density_f /= 1.0 - de_vol_f;
		}
	}
	
	return 0;
}
