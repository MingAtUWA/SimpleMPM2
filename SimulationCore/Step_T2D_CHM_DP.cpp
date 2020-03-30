#include "SimulationCore_pcp.h"

#include <cmath>
#include "Step_T2D_CHM_DP.h"

Step_T2D_CHM_DP::Step_T2D_CHM_DP() :
	Step(&solve_substep_T2D_CHM_DP), model(nullptr) {}

Step_T2D_CHM_DP::~Step_T2D_CHM_DP() {}

int Step_T2D_CHM_DP::init_calculation(void)
{
	Model_T2D_CHM_DP &md = *model;

	for (size_t pcl_id = 0; pcl_id < md.spcl_num; ++pcl_id)
	{
		SolidParticle &pcl = md.spcls[pcl_id];
		pcl.pe = (Element *)1;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
	}

	for (size_t pcl_id = 0; pcl_id < md.fpcl_num; ++pcl_id)
	{
		FluidParticle &pcl = md.fpcls[pcl_id];
		pcl.pe = (Element *)1;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
	}

	return 0;
}

int Step_T2D_CHM_DP::finalize_calculation(void) { return 0; }

namespace
{
// helper functions
bool fluid_pcl_is_in_solid(Model_T2D_CHM_DP::FluidParticle &fpcl)
{
	return fpcl.in_solid && fpcl.pe->spcls;
}

}

int solve_substep_T2D_CHM_DP(void *_self)
{
	typedef Model_T2D_CHM_DP::Node Node;
	typedef Model_T2D_CHM_DP::Element Element;
	typedef Model_T2D_CHM_DP::SolidParticle SolidParticle;
	typedef Model_T2D_CHM_DP::FluidParticle FluidParticle;
	Step_T2D_CHM_DP &self = *(Step_T2D_CHM_DP *)(_self);
	Model_T2D_CHM_DP &md = *self.model;

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		// solid
		n.has_mp_s = false;
		n.m_s = 0.0;
		n.vx_s = 0.0;
		n.vy_s = 0.0;
		n.fx_ext_s = 0.0;
		n.fy_ext_s = 0.0;
		n.fx_int_s = 0.0;
		n.fy_int_s = 0.0;
		// fluid phase
		n.has_mp_f = false;
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
	}

	// init elements
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		e.spcls = nullptr;
		e.fpcls = nullptr;
		e.vol_ps = 0.0;
		e.vol_pf = 0.0;
		e.vol_m  = 0.0;
		e.n = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
		e.p = 0.0;
	}

	// init particles
	// detect fluid in solid ...
	for (size_t pcl_id = 0; pcl_id < md.spcl_num; ++pcl_id)
	{
		SolidParticle &pcl = md.spcls[pcl_id];
		if (pcl.pe)
		{
			pcl.pe = md.find_in_which_element(pcl);
			if (!pcl.pe) continue;
			pcl.pe->add_pcl(pcl);

			double vol_gs = pcl.m / pcl.density;
			pcl.vol = vol_gs / (1.0 - pcl.n);
			Element &e = *pcl.pe;
			e.vol_ps += pcl.vol;
			e.n += vol_gs;
			e.s11 += pcl.vol * pcl.s11;
			e.s22 += pcl.vol * pcl.s22;
			e.s12 += pcl.vol * pcl.s12;

			double mvx = pcl.m * pcl.vx;
			double mvy = pcl.m * pcl.vy;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.has_mp_s = true;
			n1.m_s  += pcl.N1 * pcl.m;
			n1.vx_s += pcl.N1 * mvx;
			n1.vy_s += pcl.N1 * mvy;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.has_mp_s = true;
			n2.m_s  += pcl.N2 * pcl.m;
			n2.vx_s += pcl.N2 * mvx;
			n2.vy_s += pcl.N2 * mvy;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.has_mp_s = true;
			n3.m_s  += pcl.N3 * pcl.m;
			n3.vx_s += pcl.N3 * mvx;
			n3.vy_s += pcl.N3 * mvy;
		}
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		// void of solid phase
		double vol_sv = e.vol_ps - e.n;
		e.n = vol_sv / e.vol_ps;
	}

	for (size_t pcl_id = 0; pcl_id < md.fpcl_num; ++pcl_id)
	{
		FluidParticle &pcl = md.fpcls[pcl_id];
		if (pcl.pe)
		{
			pcl.pe = md.find_in_which_element(pcl);
			if (!pcl.pe) continue;
			pcl.pe->add_pcl(pcl);

			Element &e = *pcl.pe;
			double vol_f = pcl.m / pcl.density;
			if (fluid_pcl_is_in_solid(pcl))
			{
				pcl.vol = vol_f / e.n;
				e.vol_m += pcl.vol;
			}
			else
			{
				pcl.vol = vol_f;
				e.vol_pf += pcl.vol;
			}
			e.p += pcl.vol * pcl.p;

			// map velocity
			double mvx = pcl.m * pcl.vx;
			double mvy = pcl.m * pcl.vy;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.m_f += pcl.N1 * pcl.m;
			n1.vx_f += pcl.N1 * mvx;
			n1.vy_f += pcl.N1 * mvy;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.m_f += pcl.N2 * pcl.m;
			n2.vx_f += pcl.N2 * mvx;
			n2.vy_f += pcl.N2 * mvy;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.m_f += pcl.N3 * pcl.m;
			n3.vx_f += pcl.N3 * mvx;
			n3.vy_f += pcl.N3 * mvy;
		}
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.vol_m > e.vol_ps)
		{
			e.vol_pf += e.vol_m - e.vol_ps;
			e.vol_m = e.vol_ps;
		}
	}

	// update nodal momentum
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp_s)
		{
			n.vx_s /= n.m_s;
			n.vy_s /= n.m_s;
		}
		if (n.has_mp_f)
		{
			n.vx_f /= n.m_f;
			n.vy_f /= n.m_f;
		}
	}

	// apply velocity bc
	for (size_t v_id = 0; v_id < md.svx_num; ++v_id)
	{
		Node &n = md.nodes[md.svxs[v_id].node_id];
		n.vx_s = md.svxs[v_id].v;
	}
	for (size_t v_id = 0; v_id < md.svy_num; ++v_id)
	{
		Node &n = md.nodes[md.svys[v_id].node_id];
		n.vy_s = md.svys[v_id].v;
	}
	for (size_t v_id = 0; v_id < md.fvx_num; ++v_id)
	{
		Node &n = md.nodes[md.fvxs[v_id].node_id];
		n.vx_f = md.fvxs[v_id].v;
	}
	for (size_t v_id = 0; v_id < md.fvy_num; ++v_id)
	{
		Node &n = md.nodes[md.fvys[v_id].node_id];
		n.vy_f = md.fvys[v_id].v;
	}

	//double de_vol_f_rate;
	//for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	//{
	//	Element &e = md.elems[e_id];
	//	if (e.spcls)
	//	{
	//		e.s11 /= e.vol;
	//		e.s22 /= e.vol;
	//		e.s12 /= e.vol;
	//	}

	//	if (e.fpcls)
	//	{
	//		e.p /= e.vol_pf;
	//	}

	//	Node &n1 = md.nodes[e.n1];
	//	Node &n2 = md.nodes[e.n2];
	//	Node &n3 = md.nodes[e.n3];

	//	// volume strain rate
	//	de_vol_f_rate = -(1.0 - e.n) / e.n *
	//					(n1.vx_s * e.dN1_dx + n2.vx_s * e.dN2_dx + n3.vx_s * e.dN3_dx
	//					+ n1.vy_s * e.dN1_dy + n2.vy_s * e.dN2_dy + n3.vy_s * e.dN3_dy)
	//					-(n1.vx_f * e.dN1_dx + n2.vx_f * e.dN2_dx + n3.vx_f * e.dN3_dx
	//					+ n1.vy_f * e.dN1_dy + n2.vy_f * e.dN2_dy + n3.vy_f * e.dN3_dy);
	//	// add bulk viscosity
	//	e.p -= self.bv_ratio * de_vol_f_rate;
	//		
	//	// node 1
	//	n1.fx_int_s += (e.dN1_dx * (e.s11 - (1.0 - e.n) * e.p) + e.dN1_dy * e.s12) * e.vol;
	//	n1.fy_int_s += (e.dN1_dx * e.s12 + e.dN1_dy * (e.s22 - (1.0 - e.n) * e.p)) * e.vol;
	//	n1.fx_int_f += (e.dN1_dx * e.n * -e.p) * e.vol;
	//	n1.fy_int_f += (e.dN1_dy * e.n * -e.p) * e.vol;
	//	// node 2
	//	n2.fx_int_s += (e.dN2_dx * (e.s11 - (1.0 - e.n) * e.p) + e.dN2_dy * e.s12) * e.vol;
	//	n2.fy_int_s += (e.dN2_dx * e.s12 + e.dN2_dy * (e.s22 - (1.0 - e.n) * e.p)) * e.vol;
	//	n2.fx_int_f += (e.dN2_dx * e.n * -e.p) * e.vol;
	//	n2.fy_int_f += (e.dN2_dy * e.n * -e.p) * e.vol;
	//	// node 3
	//	n3.fx_int_s += (e.dN3_dx * (e.s11 - (1.0 - e.n) * e.p) + e.dN3_dy * e.s12) * e.vol;
	//	n3.fy_int_s += (e.dN3_dx * e.s12 + e.dN3_dy * (e.s22 - (1.0 - e.n) * e.p)) * e.vol;
	//	n3.fx_int_f += (e.dN3_dx * e.n * -e.p) * e.vol;
	//	n3.fy_int_f += (e.dN3_dy * e.n * -e.p) * e.vol;
	//}

	// seepage force
	for (size_t pcl_id = 0; pcl_id < md.fpcl_num; ++pcl_id)
	{
		FluidParticle &pcl = md.fpcls[pcl_id];


	}

	// body force
	double bf_mag;
	for (size_t bf_id = 0; bf_id < md.sbfx_num; ++bf_id)
	{
		BodyForce &bf = md.sbfxs[bf_id];
		SolidParticle &pcl = md.spcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			bf_mag = pcl.m * bf.bf;
			Node &n1 = md.nodes[e.n1];
			n1.fx_ext_s += pcl.N1 * bf_mag;
			Node &n2 = md.nodes[e.n2];
			n2.fx_ext_s += pcl.N2 * bf_mag;
			Node &n3 = md.nodes[e.n3];
			n3.fx_ext_s += pcl.N3 * bf_mag;
		}
	}
	for (size_t bf_id = 0; bf_id < md.sbfy_num; ++bf_id)
	{
		BodyForce &bf = md.sbfys[bf_id];
		SolidParticle &pcl = md.spcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			bf_mag = pcl.m * bf.bf;
			Node &n1 = md.nodes[e.n1];
			n1.fy_ext_s += pcl.N1 * bf_mag;
			Node &n2 = md.nodes[e.n2];
			n2.fy_ext_s += pcl.N2 * bf_mag;
			Node &n3 = md.nodes[e.n3];
			n3.fy_ext_s += pcl.N3 * bf_mag;
		}
	}
	for (size_t bf_id = 0; bf_id < md.fbfx_num; ++bf_id)
	{
		BodyForce &bf = md.fbfxs[bf_id];
		FluidParticle &pcl = md.fpcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			bf_mag = pcl.m * bf.bf;
			Node &n1 = md.nodes[e.n1];
			n1.fx_ext_f += pcl.N1 * bf_mag;
			Node &n2 = md.nodes[e.n2];
			n2.fx_ext_f += pcl.N2 * bf_mag;
			Node &n3 = md.nodes[e.n3];
			n3.fx_ext_f += pcl.N3 * bf_mag;
		}
	}
	for (size_t bf_id = 0; bf_id < md.fbfy_num; ++bf_id)
	{
		BodyForce &bf = md.fbfys[bf_id];
		FluidParticle &pcl = md.fpcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			bf_mag = pcl.m * bf.bf;
			Node &n1 = md.nodes[e.n1];
			n1.fy_ext_f += pcl.N1 * bf_mag;
			Node &n2 = md.nodes[e.n2];
			n2.fy_ext_f += pcl.N2 * bf_mag;
			Node &n3 = md.nodes[e.n3];
			n3.fy_ext_f += pcl.N3 * bf_mag;
		}
	}

	// surface force
	for (size_t tf_id = 0; tf_id < md.tx_num; ++tf_id)
	{
		TractionBC_MPM &tf = md.txs[tf_id];
		SolidParticle &pcl = md.spcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			n1.fx_ext_s += pcl.N1 * tf.t;
			Node &n2 = md.nodes[e.n2];
			n2.fx_ext_s += pcl.N2 * tf.t;
			Node &n3 = md.nodes[e.n3];
			n3.fx_ext_s += pcl.N3 * tf.t;
		}
	}
	for (size_t tf_id = 0; tf_id < md.ty_num; ++tf_id)
	{
		TractionBC_MPM &tf = md.tys[tf_id];
		SolidParticle &pcl = md.spcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			n1.fy_ext_s += pcl.N1 * tf.t;
			Node &n2 = md.nodes[e.n2];
			n2.fy_ext_s += pcl.N2 * tf.t;
			Node &n3 = md.nodes[e.n3];
			n3.fy_ext_s += pcl.N3 * tf.t;
		}
	}

	// update nodal acceleration of fluid pahse
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp_s)
		{
			n.ax_s = n.fx_ext_s - n.fx_int_s + n.fx_drag / n.m_s;
			n.ay_s = n.fy_ext_s - n.fy_int_s + n.fy_drag / n.m_s;
		}

		if (n.has_mp_f)
		{
			n.ax_f = n.fx_ext_f - n.fx_int_f - n.fx_drag / n.m_f;
			n.ay_f = n.fy_ext_f - n.fy_int_f - n.fy_drag / n.m_f;
		}
	}

	// apply acceleration bc
	for (size_t a_id = 0; a_id < md.sax_num; ++a_id)
	{
		Node &n = md.nodes[md.saxs[a_id].node_id];
		n.ax_s = md.saxs[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.say_num; ++a_id)
	{
		Node &n = md.nodes[md.says[a_id].node_id];
		n.ay_s = md.says[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.fax_num; ++a_id)
	{
		Node &n = md.nodes[md.faxs[a_id].node_id];
		n.ax_f = md.faxs[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.fay_num; ++a_id)
	{
		Node &n = md.nodes[md.fays[a_id].node_id];
		n.ay_f = md.fays[a_id].a;
	}

	// update nodal momentum
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp_s)
		{
			n.vx_s += n.ax_s * self.dtime;
			n.vy_s += n.ay_s * self.dtime;
		}

		if (n.has_mp_f)
		{
			n.vx_f += n.ax_f * self.dtime;
			n.vy_f += n.ay_f * self.dtime;
		}
	}

	// apply velocity bc
	for (size_t v_id = 0; v_id < md.svx_num; ++v_id)
	{
		Node &n = md.nodes[md.svxs[v_id].node_id];
		n.vx_s = md.svxs[v_id].v;
		n.ax_s = 0.0;
	}
	for (size_t v_id = 0; v_id < md.svy_num; ++v_id)
	{
		Node &n = md.nodes[md.svys[v_id].node_id];
		n.vy_s = md.svys[v_id].v;
		n.ay_s = 0.0;
	}
	for (size_t v_id = 0; v_id < md.fvx_num; ++v_id)
	{
		Node &n =  md.nodes[md.fvxs[v_id].node_id];
		n.vx_f = md.fvxs[v_id].v;
		n.ax_f = 0.0;
	}
	for (size_t v_id = 0; v_id < md.fvy_num; ++v_id)
	{
		Node &n = md.nodes[md.fvys[v_id].node_id];
		n.vy_f = md.fvys[v_id].v;
		n.ay_f = 0.0;
	}

	// update displacement increment of both phases
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp_s)
		{
			n.dux_s = n.vx_s * self.dtime;
			n.duy_s = n.vy_s * self.dtime;
		}
		if (n.has_mp_f)
		{
			n.dux_f = n.vx_f * self.dtime;
			n.duy_f = n.vy_f * self.dtime;
		}
	}

	// map variables back to particles and update their variables
	double de11, de22, de12, de_vol_s;
	double de_vol_pf, de_vol_m;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		Node &n1 = md.nodes[e.n1];
		Node &n2 = md.nodes[e.n2];
		Node &n3 = md.nodes[e.n3];

		if (e.spcls)
		{
			// strain increment
			de11 = n1.dux_s * e.dN1_dx + n2.dux_s * e.dN2_dx + n3.dux_s * e.dN3_dx;
			de22 = n1.duy_s * e.dN1_dy + n2.duy_s * e.dN2_dy + n3.duy_s * e.dN3_dy;
			de12 = (n1.dux_s * e.dN1_dy + n2.dux_s * e.dN2_dy + n3.dux_s * e.dN3_dy
				+ n1.duy_s * e.dN1_dx + n2.duy_s * e.dN2_dx + n3.duy_s * e.dN3_dx) * 0.5;

			// volumetric strain of solid phase
			de_vol_s = de11 + de22;
		}

		if (e.fpcls)
		{
			//de_vol_pf = ;
			de_vol_m = (1.0 - e.n) / e.n * -de_vol_s
				- (n1.dux_f * e.dN1_dx + n2.dux_f * e.dN2_dx + n3.dux_f * e.dN3_dx
					+ n1.duy_f * e.dN1_dy + n2.duy_f * e.dN2_dy + n3.duy_f * e.dN3_dy);
		}
	}

	for (size_t pcl_id = 0; pcl_id < md.spcl_num; ++pcl_id)
	{
		//SolidParticle &pcl = md.spcls[pcl_id];

		//// velocity
		//pcl.vx_s += (n1.ax_s * pcl.N1 + n2.ax_s * pcl.N2 + n3.ax_s * pcl.N3) * self.dtime;
		//pcl.vy_s += (n1.ay_s * pcl.N1 + n2.ay_s * pcl.N2 + n3.ay_s * pcl.N3) * self.dtime;
		//pcl.vx_f += (n1.ax_f * pcl.N1 + n2.ax_f * pcl.N2 + n3.ax_f * pcl.N3) * self.dtime;
		//pcl.vy_f += (n1.ay_f * pcl.N1 + n2.ay_f * pcl.N2 + n3.ay_f * pcl.N3) * self.dtime;

		//// displacement
		//pcl.ux_s += n1.dux_s * pcl.N1 + n2.dux_s * pcl.N2 + n3.dux_s * pcl.N3;
		//pcl.uy_s += n1.duy_s * pcl.N1 + n2.duy_s * pcl.N2 + n3.duy_s * pcl.N3;
		//pcl.ux_f += n1.dux_f * pcl.N1 + n2.dux_f * pcl.N2 + n3.dux_f * pcl.N3;
		//pcl.uy_f += n1.duy_f * pcl.N1 + n2.duy_f * pcl.N2 + n3.duy_f * pcl.N3;

		//// update position
		//pcl.x = pcl.x_ori + pcl.ux_s;
		//pcl.y = pcl.y_ori + pcl.uy_s;

		//// strain
		//pcl.e11 += de11;
		//pcl.e22 += de22;
		//pcl.e12 += de12;

		//// stress
		//// directly assign element-wise stress
		//pcl.s11 = e.s11;
		//pcl.s22 = e.s22;
		//pcl.s12 = e.s12;
		//pcl.p = e.p;

		//// fluid density
		//pcl.density_f /= 1.0 - de_vol_f;

		//// porosity
		//pcl.n = (de_vol_s + pcl.n) / (1.0 + de_vol_s);
	}
	
	for (size_t pcl_id = 0; pcl_id < md.fpcl_num; ++pcl_id)
	{
		FluidParticle &pcl = md.fpcls[pcl_id];
		

	}

	return 0;
}
