#include "SimulationCore_pcp.h"

#include "Step_S2D_ME_s_Geostatic.h"

//size_t Step_S2D_ME_s_Geostatic::distribute_pcl_mat_to_elems(Particle &pcl)
//{
//	Model_S2D_ME_s &md = *static_cast<Model_S2D_ME_s *>(model);
//
//	// half length
//	pcl.vol = pcl.m / pcl.density;
//	double hlen = sqrt(pcl.vol/pcl.ar) * 0.5;
//	double wlen = hlen * pcl.ar;
//	double xl = pcl.x - wlen;
//	double xu = pcl.x + wlen;
//	double yl = pcl.y - hlen;
//	double yu = pcl.y + hlen;
//
//	// check if particles is out of mesh
//	if (xu <= md.x0 || xl >= md.xn || yu <= md.y0 || yl >= md.yn)
//		return 0;
//
//	// trim particle if it lies at edge
//	if (xl < md.x0)
//		xl = md.x0;
//	if (xu > md.xn)
//		xu = md.xn;
//	if (yl < md.y0)
//		yl = md.y0;
//	if (yu > md.yn)
//		yu = md.yn;
//
//	size_t xl_id = size_t(floor((xl - md.x0 + md.h_tol) / md.hx));
//	size_t xu_id = size_t(ceil((xu - md.x0 - md.h_tol) / md.hx));
//	size_t x_num = xu_id - xl_id;
//	size_t yl_id = size_t(floor((yl - md.y0 + md.h_tol) / md.hy));
//	size_t yu_id = size_t(ceil((yu - md.y0 - md.h_tol) / md.hy));
//	size_t y_num = yu_id - yl_id;
//	size_t elem_num = x_num * y_num;
//
//	Element *pe;
//	if (elem_num == 1)
//	{
//		pe = md.elems + md.elem_x_num * yl_id + xl_id;
//		pe->ve_vol += (xu - xl) * (yu - yl);
//		return 1;
//	}
//
//	double x_len1, x_len2, y_len1, y_len2;
//	if (x_num == 1)
//	{
//		if (y_num == 2)
//		{
//			x_len1 = xu - xl;
//			y_len1 = md.y0 + double(yl_id + 1) * md.hy - yl;
//			y_len2 = yu - yl - y_len1;
//			// elem 1
//			pe = md.elems + md.elem_x_num * yl_id + xl_id;
//			pe->ve_vol += x_len1 * y_len1;
//			// elem 2
//			pe += md.elem_x_num;
//			pe->ve_vol += x_len1 * y_len2;
//			return 2;
//		}
//	}
//	else if (x_num == 2)
//	{
//		x_len1 = md.x0 + double(xl_id + 1) * md.hx - xl;
//		x_len2 = xu - xl - x_len1;
//		if (y_num == 1)
//		{
//			y_len1 = yu - yl;
//			// elem 1
//			pe = md.elems + md.elem_x_num * yl_id + xl_id;
//			pe->ve_vol += x_len1 * y_len1;
//			// elem 2
//			++pe;
//			pe->ve_vol += x_len2 * y_len1;
//			return 2;
//		}
//		else if (y_num == 2)
//		{
//			y_len1 = md.y0 + double(yl_id + 1) * md.hy - yl;
//			y_len2 = yu - yl - y_len1;
//			// elem 1
//			pe = md.elems + md.elem_x_num * yl_id + xl_id;
//			pe->ve_vol += x_len1 * y_len1;
//			// elem 2
//			++pe;
//			pe->ve_vol += x_len2 * y_len1;
//			// elem 4
//			pe += md.elem_x_num;
//			pe->ve_vol += x_len2 * y_len2;
//			// elem 3
//			--pe;
//			pe->ve_vol += x_len1 * y_len2;
//			return 4;
//		}
//	}
//
//	double *pcl_x_len, *pcl_y_len;
//	// x
//	pcl_x_len_buf.reset();
//	if (x_num == 1)
//	{
//		pcl_x_len = &x_len1;
//		x_len1 = xu - xl;
//	}
//	else
//	{
//		pcl_x_len = pcl_x_len_buf.alloc(x_num);
//		pcl_x_len[0] = md.x0 + double(xl_id + 1) * md.hx - xl;
//		for (size_t i = 1; i < x_num-1; ++i)
//			pcl_x_len[i] = md.hx;
//		pcl_x_len[x_num-1] = xu - double(xu_id - 1) * md.hx - md.x0;
//	}
//	// y
//	pcl_y_len_buf.reset();
//	if (y_num == 1)
//	{
//		pcl_y_len = &y_len1;
//		y_len1 = yu - yl;
//	}
//	else
//	{
//		pcl_y_len = pcl_y_len_buf.alloc(y_num);
//		pcl_y_len[0] = md.y0 + double(yl_id + 1) * md.hy - yl;
//		for (size_t i = 1; i < y_num-1; ++i)
//			pcl_y_len[i] = md.hy;
//		pcl_y_len[y_num-1] = yu - double(yu_id - 1) * md.hy - md.y0;
//	}
//
//	pe = md.elems + md.elem_x_num * yl_id + xl_id;
//	for (size_t y_id = 0; y_id < y_num; ++y_id)
//	{
//		for (size_t x_id = 0; x_id < x_num; ++x_id)
//			pe[x_id].ve_vol += pcl_x_len[x_id] * pcl_y_len[y_id];
//		pe += md.elem_x_num;
//	}
//
//	return elem_num;
//}

int solve_substep_S2D_ME_s_Geostatic_VE(void *_self)
{
	typedef Step_S2D_ME_s_Geostatic::ShapeFuncValue ShapeFuncValue;
	typedef Step_S2D_ME_s_Geostatic::Node Node;
	typedef Step_S2D_ME_s_Geostatic::Element Element;
	typedef Step_S2D_ME_s_Geostatic::Particle Particle;
	
	Step_S2D_ME_s_Geostatic &self = *static_cast<Step_S2D_ME_s_Geostatic *>(_self);
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
		// mixed integration
		e.mi_vol = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
		// volume enhancement
		e.ve_vol = 0.0;
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

			// volume enhancement
			md.distribute_pcl_mat_to_elems(pcl);

			// map velocity
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
			// average stress
			e.s11 /= e.mi_vol;
			e.s22 /= e.mi_vol;
			e.s12 /= e.mi_vol;
			// internal force
			// node 1
			Node &n1 = *e.pn1;
			n1.fx_int += (elem_sfv.dN1_dx * e.s11 + elem_sfv.dN1_dy * e.s12) * e.ve_vol;
			n1.fy_int += (elem_sfv.dN1_dx * e.s12 + elem_sfv.dN1_dy * e.s22) * e.ve_vol;
			// node 2
			Node &n2 = *e.pn2;
			n2.fx_int += (elem_sfv.dN2_dx * e.s11 + elem_sfv.dN2_dy * e.s12) * e.ve_vol;
			n2.fy_int += (elem_sfv.dN2_dx * e.s12 + elem_sfv.dN2_dy * e.s22) * e.ve_vol;
			// node 3
			Node &n3 = *e.pn3;
			n3.fx_int += (elem_sfv.dN3_dx * e.s11 + elem_sfv.dN3_dy * e.s12) * e.ve_vol;
			n3.fy_int += (elem_sfv.dN3_dx * e.s12 + elem_sfv.dN3_dy * e.s22) * e.ve_vol;
			// node 4
			Node &n4 = *e.pn4;
			n4.fx_int += (elem_sfv.dN4_dx * e.s11 + elem_sfv.dN4_dy * e.s12) * e.ve_vol;
			n4.fy_int += (elem_sfv.dN4_dx * e.s12 + elem_sfv.dN4_dy * e.s22) * e.ve_vol;
		}
	}

	//for (size_t y_id = 0; y_id < md.elem_y_num; ++y_id)
	//{
	//	for (size_t x_id = 0; x_id < md.elem_x_num; ++x_id)
	//	{
	//		Element &e = md.elems[y_id * md.elem_x_num + x_id];
	//		std::cout << e.ve_vol << " ";
	//	}
	//	std::cout << "\n";
	//}

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

	// update displacement increment of both phases
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.dux = n.vx * self.dtime;
			n.duy = n.vy * self.dtime;
		}
	}

	double e_kin = 0.0;
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

			// velocity
			pcl.vx += (n1.ax * sfv.N1 + n2.ax * sfv.N2
				+ n3.ax * sfv.N3 + n4.ax * sfv.N4) * self.dtime;
			pcl.vy += (n1.ay * sfv.N1 + n2.ay * sfv.N2
				+ n3.ay * sfv.N3 + n4.ay * sfv.N4) * self.dtime;

			e_kin += pcl.m * (pcl.vx * pcl.vx + pcl.vy * pcl.vy);
		}
	}

	if (e_kin < self.prev_e_kin)
	{
		for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
		{
			Particle &pcl = md.pcls[pcl_id];
			pcl.vx = 0.0;
			pcl.vy = 0.0;
		}
		self.prev_e_kin = 0.0;
		return 0;
	}
	self.prev_e_kin = e_kin;

	double de11, de22, de12;
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
		}
	}

	return 0;
}
