#include "SimulationCore_pcp.h"

#include "Step_S2D_ME_s.h"

Step_S2D_ME_s::Step_S2D_ME_s() :
	Step(&solve_substep_S2D_ME_s) {}

Step_S2D_ME_s::~Step_S2D_ME_s() {}

int Step_S2D_ME_s::init_calculation(void)
{
	Model_S2D_ME_s &md = *static_cast<Model_S2D_ME_s *>(model);

	if (is_first_step) {}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		pcl.pe = md.elems;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
	}

	return 0;
}

int Step_S2D_ME_s::finalize_calculation(void)
{
	Model_S2D_ME_s &md = *static_cast<Model_S2D_ME_s *>(model);
	
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		Particle &pcl = md.pcls[p_id];
		pcl.vx = 0.0;
		pcl.vy = 0.0;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
	}

	return 0;
}

int solve_substep_S2D_ME_s(void *_self)
{
	typedef Step_S2D_ME_s::ShapeFuncValue ShapeFuncValue;
	typedef Step_S2D_ME_s::Node Node;
	typedef Step_S2D_ME_s::Element Element;
	typedef Step_S2D_ME_s::Particle Particle;
	
	Step_S2D_ME_s &self = *static_cast<Step_S2D_ME_s *>(_self);
	Model_S2D_ME_s &md = *static_cast<Model_S2D_ME_s *>(self.model);
	ShapeFuncValue &elem_sfv = md.elem_sfv;

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		n.has_mp = false;
		n.m = 0.0;
		n.vx = 0.0;
		n.vy = 0.0;
		n.fx_ext = 0.0;
		n.fy_ext = 0.0;
		n.fx_int = 0.0;
		n.fy_int = 0.0;
	}

	// init elements
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		e.has_pcl = false;
		e.mi_vol = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
	}
	
	// init particles
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			if (!md.find_in_which_element(pcl))
				continue;

			pcl.vol = pcl.m / pcl.density;
			
			// mixed integration
			Element &e = *pcl.pe;
			e.has_pcl = true;
			e.mi_vol += pcl.vol;
			e.s11 += pcl.s11 * pcl.vol;
			e.s22 += pcl.s22 * pcl.vol;
			e.s12 += pcl.s12 * pcl.vol;

			ShapeFuncValue &sfv = pcl.sfv;
			double mvx = pcl.m * pcl.vx;
			double mvy = pcl.m * pcl.vy;
			// node 1
			Node &n1 = *e.pn1;
			n1.has_mp = true;
			n1.m += sfv.N1 * pcl.m;
			n1.vx += sfv.N1 * mvx;
			n1.vy += sfv.N1 * mvy;
			// node 2
			Node &n2 = *e.pn2;
			n2.has_mp = true;
			n2.m += sfv.N2 * pcl.m;
			n2.vx += sfv.N2 * mvx;
			n2.vy += sfv.N2 * mvy;
			// node 3
			Node &n3 = *e.pn3;
			n3.has_mp = true;
			n3.m += sfv.N3 * pcl.m;
			n3.vx += sfv.N3 * mvx;
			n3.vy += sfv.N3 * mvy;
			// node 4
			Node &n4 = *e.pn4;
			n4.has_mp = true;
			n4.m += sfv.N4 * pcl.m;
			n4.vx += sfv.N4 * mvx;
			n4.vy += sfv.N4 * mvy;
		}
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.has_pcl)
		{
			e.s11 /= e.mi_vol;
			e.s22 /= e.mi_vol;
			e.s12 /= e.mi_vol;
		}
		if (e.mi_vol > md.elem_area)
			e.mi_vol = md.elem_area;

		Node &n1 = *e.pn1;
		n1.fx_int += (elem_sfv.dN1_dx * e.s11 + elem_sfv.dN1_dy * e.s12) * e.mi_vol;
		n1.fy_int += (elem_sfv.dN1_dx * e.s12 + elem_sfv.dN1_dy * e.s22) * e.mi_vol;
		Node &n2 = *e.pn2;
		n2.fx_int += (elem_sfv.dN2_dx * e.s11 + elem_sfv.dN2_dy * e.s12) * e.mi_vol;
		n2.fy_int += (elem_sfv.dN2_dx * e.s12 + elem_sfv.dN2_dy * e.s22) * e.mi_vol;
		Node &n3 = *e.pn3;
		n3.fx_int += (elem_sfv.dN3_dx * e.s11 + elem_sfv.dN3_dy * e.s12) * e.mi_vol;
		n3.fy_int += (elem_sfv.dN3_dx * e.s12 + elem_sfv.dN3_dy * e.s22) * e.mi_vol;
		Node &n4 = *e.pn4;
		n4.fx_int += (elem_sfv.dN4_dx * e.s11 + elem_sfv.dN4_dy * e.s12) * e.mi_vol;
		n4.fy_int += (elem_sfv.dN4_dx * e.s12 + elem_sfv.dN4_dy * e.s22) * e.mi_vol;
	}

	// body force
	for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
	{
		BodyForce &bf = md.bfxs[bf_id];
		Particle &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			ShapeFuncValue &sfv = pcl.sfv;
			// node 1
			Node &n1 = *e.pn1;
			n1.fx_ext += sfv.N1 * bf.bf;
			// node 2
			Node &n2 = *e.pn2;
			n2.fx_ext += sfv.N2 * bf.bf;
			// node 3
			Node &n3 = *e.pn3;
			n3.fx_ext += sfv.N3 * bf.bf;
			// node 4
			Node &n4 = *e.pn4;
			n4.fx_ext += sfv.N4 * bf.bf;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
	{
		BodyForce &bf = md.bfys[bf_id];
		Particle &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			ShapeFuncValue &sfv = pcl.sfv;
			// node 1
			Node &n1 = *e.pn1;
			n1.fy_ext += sfv.N1 * bf.bf;
			// node 2
			Node &n2 = *e.pn2;
			n2.fy_ext += sfv.N2 * bf.bf;
			// node 3
			Node &n3 = *e.pn3;
			n3.fy_ext += sfv.N3 * bf.bf;
			// node 4
			Node &n4 = *e.pn4;
			n4.fy_ext += sfv.N4 * bf.bf;
		}
	}

	// surface force
	for (size_t tf_id = 0; tf_id < md.tx_num; ++tf_id)
	{
		TractionBC_MPM &tf = md.txs[tf_id];
		Particle &pcl = md.pcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			ShapeFuncValue &sfv = pcl.sfv;
			// node 1
			Node &n1 = *e.pn1;
			n1.fx_ext += sfv.N1 * tf.t;
			// node 2
			Node &n2 = *e.pn2;
			n2.fx_ext += sfv.N2 * tf.t;
			// node 3
			Node &n3 = *e.pn3;
			n3.fx_ext += sfv.N3 * tf.t;
			// node 4
			Node &n4 = *e.pn4;
			n4.fx_ext += sfv.N4 * tf.t;
		}
	}
	for (size_t tf_id = 0; tf_id < md.ty_num; ++tf_id)
	{
		TractionBC_MPM &tf = md.tys[tf_id];
		Particle &pcl = md.pcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			ShapeFuncValue &sfv = pcl.sfv;
			// node 1
			Node &n1 = *e.pn1;
			n1.fy_ext += sfv.N1 * tf.t;
			// node 2
			Node &n2 = *e.pn2;
			n2.fy_ext += sfv.N2 * tf.t;
			// node 3
			Node &n3 = *e.pn3;
			n3.fy_ext += sfv.N3 * tf.t;
			// node 4
			Node &n4 = *e.pn4;
			n4.fy_ext += sfv.N4 * tf.t;
		}
	}

	// update nodal acceleration of fluid pahse
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.ax = (n.fx_ext - n.fx_int) / n.m;
			n.ay = (n.fy_ext - n.fy_int) / n.m;
		}
	}
	for (size_t a_id = 0; a_id < md.ax_num; ++a_id)
	{
		AccelerationBC &abc = md.axs[a_id];
		md.nodes[abc.node_id].ax = abc.a;
	}
	for (size_t a_id = 0; a_id < md.ay_num; ++a_id)
	{
		AccelerationBC &abc = md.ays[a_id];
		md.nodes[abc.node_id].ay = abc.a;
	}

	// update nodal momentum of fluid phase
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.vx /= n.m;
			n.vx += n.ax * self.dtime;
			n.vy /= n.m;
			n.vy += n.ay * self.dtime;
		}
	}
	// apply velocity boundary conditions of fluid phase
	for (size_t v_id = 0; v_id < md.vx_num; ++v_id)
	{
		VelocityBC &vbc = md.vxs[v_id];
		Node &n = md.nodes[vbc.node_id];
		n.vx = vbc.v;
		n.ax = 0.0;
	}
	for (size_t v_id = 0; v_id < md.vy_num; ++v_id)
	{
		VelocityBC &vbc = md.vys[v_id];
		Node &n = md.nodes[vbc.node_id];
		n.vy = vbc.v;
		n.ay = 0.0;
	}

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.dux = n.vx * self.dtime;
			n.duy = n.vy * self.dtime;
		}
	}

	double de11, de22, de12, de_vol;
	double ds11, ds22, ds12;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = *e.pn1;
			Node &n2 = *e.pn2;
			Node &n3 = *e.pn3;
			Node &n4 = *e.pn4;
			ShapeFuncValue &sfv = pcl.sfv;

			// position
			pcl.ux += n1.dux * sfv.N1 + n2.dux * sfv.N2
					+ n3.dux * sfv.N3 + n4.dux * sfv.N4;
			pcl.uy += n1.duy * sfv.N1 + n2.duy * sfv.N2
					+ n3.duy * sfv.N3 + n4.duy * sfv.N4;
			pcl.x = pcl.x_ori + pcl.ux;
			pcl.y = pcl.y_ori + pcl.uy;

			// velocity
			pcl.vx += (n1.ax * sfv.N1 + n2.ax * sfv.N2
					 + n3.ax * sfv.N3 + n4.ax * sfv.N4) * self.dtime;
			pcl.vy += (n1.ay * sfv.N1 + n2.ay * sfv.N2
					 + n3.ay * sfv.N3 + n4.ay * sfv.N4) * self.dtime;

			// strain increment
			de11 = n1.dux * sfv.dN1_dx + n2.dux * sfv.dN2_dx
				 + n3.dux * sfv.dN3_dx + n4.dux * sfv.dN4_dx;
			de22 = n1.duy * sfv.dN1_dy + n2.duy * sfv.dN2_dy
				 + n3.duy * sfv.dN3_dy + n4.duy * sfv.dN4_dy;
			de12 = (n1.dux * sfv.dN1_dy + n2.dux * sfv.dN2_dy
				  + n3.dux * sfv.dN3_dy + n4.dux * sfv.dN4_dy
				  + n1.duy * sfv.dN1_dx + n2.duy * sfv.dN2_dx
				  + n3.duy * sfv.dN3_dx + n4.duy * sfv.dN4_dx) * 0.5;
			pcl.e11 += de11;
			pcl.e22 += de22;
			pcl.e12 += de12;

			// update stress
			double E_tmp = md.E / (1.0 + md.niu) / (1.0 - 2.0 * md.niu);
			ds11 = E_tmp * ((1.0 - md.niu) * de11 + md.niu * de22);
			ds22 = E_tmp * (md.niu * de11 + (1.0 - md.niu) * de22);
			ds12 = 2.0 * md.E / (2.0 * (1.0 + md.niu)) * de12;
			pcl.s11 += ds11;
			pcl.s22 += ds22;
			pcl.s12 += ds12;

			// density
			de_vol = de11 + de22;
			pcl.density /= (1.0 + de_vol);
		}
	}

	return 0;
}
