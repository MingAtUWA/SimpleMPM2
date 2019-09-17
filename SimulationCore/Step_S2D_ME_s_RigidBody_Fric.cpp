#include "SimulationCore_pcp.h"

#include "Step_S2D_ME_s_RigidBody_Fric.h"

namespace
{
	typedef Model_S2D_ME_s_RigidBody::Particle Particle_rb;
	typedef Model_S2D_ME_s_RigidBody::Element Element_rb;
	typedef Model_S2D_ME_s_RigidBody::Node Node_rb;
};

Step_S2D_ME_s_RigidBody_Fric::Step_S2D_ME_s_RigidBody_Fric() :
	Step(&solve_substep_S2D_ME_s_RigidBody),
	h_elem_raio(0.05), h_pcl_ratio(0.1),
	model(nullptr) {}

Step_S2D_ME_s_RigidBody_Fric::~Step_S2D_ME_s_RigidBody_Fric() {}

int Step_S2D_ME_s_RigidBody_Fric::init_calculation()
{
	if (is_first_step)
	{
		for (size_t i = 0; i < model->pcl_num; ++i)
		{
			// assign non-zero value
			model->pcls[i].pelem = model->elems; 
		}
	}

	for (size_t i = 0; i < model->pcl_num; ++i)
	{
		Particle_rb &pcl = model->pcls[i];
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
	}

	model->rigid_body.init_calculation();

	return 0;
}

int Step_S2D_ME_s_RigidBody_Fric::finalize_calculation() { return 0; }

int solve_substep_S2D_ME_s_RigidBody_Fric(void *_self)
{
	Step_S2D_ME_s_RigidBody_Fric &self = *((Step_S2D_ME_s_RigidBody_Fric *)_self);
	Model_S2D_ME_s_RigidBody &md = *(self.model);

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_rb &n = md.nodes[n_id];
		n.m = 0.0;
		n.ax = 0.0;
		n.ay = 0.0;
		n.vx = 0.0;
		n.vy = 0.0;
		n.fx_ext = 0.0;
		n.fy_ext = 0.0;
		n.fx_int = 0.0;
		n.fy_int = 0.0;
		n.fx_con = 0.0;
		n.fy_con = 0.0;
	}

	// init elements
	for (size_t elem_id = 0; elem_id < md.elem_num; ++elem_id)
		md.elems[elem_id].reset();

	// map variables to nodes and cal internal force
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle_rb &pcl = md.pcls[pcl_id];
		if (pcl.pelem) // is in mesh
		{
			// find in which bg_mesh element
			if (!md.init_pcl_cal_var(pcl))
				continue;

			// add pcl to background element (for contact detection)
			pcl.pelem->add_pcl(pcl);

			// Map to nodes
			double mvx = pcl.m * pcl.vx;
			double mvy = pcl.m * pcl.vy;
			// node 1
			Node_rb &n1 = *pcl.pn1;
			n1.m  += pcl.N1 * pcl.m;
			n1.vx += pcl.N1 * mvx;
			n1.vy += pcl.N1 * mvy;
			n1.fx_int += (pcl.dN1_dx * pcl.s11 + pcl.dN1_dy * pcl.s12) * pcl.vol;
			n1.fy_int += (pcl.dN1_dx * pcl.s12 + pcl.dN1_dy * pcl.s22) * pcl.vol;
			// node 2
			Node_rb &n2 = *pcl.pn2;
			n2.m  += pcl.N2 * pcl.m;
			n2.vx += pcl.N2 * mvx;
			n2.vy += pcl.N2 * mvy;
			n2.fx_int += (pcl.dN2_dx * pcl.s11 + pcl.dN2_dy * pcl.s12) * pcl.vol;
			n2.fy_int += (pcl.dN2_dx * pcl.s12 + pcl.dN2_dy * pcl.s22) * pcl.vol;
			// node 3
			Node_rb &n3 = *pcl.pn3;
			n3.m  += pcl.N3 * pcl.m;
			n3.vx += pcl.N3 * mvx;
			n3.vy += pcl.N3 * mvy;
			n3.fx_int += (pcl.dN3_dx * pcl.s11 + pcl.dN3_dy * pcl.s12) * pcl.vol;
			n3.fy_int += (pcl.dN3_dx * pcl.s12 + pcl.dN3_dy * pcl.s22) * pcl.vol;
			// node 4
			Node_rb &n4 = *pcl.pn4;
			n4.m  += pcl.N4 * pcl.m;
			n4.vx += pcl.N4 * mvx;
			n4.vy += pcl.N4 * mvy;
			n4.fx_int += (pcl.dN4_dx * pcl.s11 + pcl.dN4_dy * pcl.s12) * pcl.vol;
			n4.fy_int += (pcl.dN4_dx * pcl.s12 + pcl.dN4_dy * pcl.s22) * pcl.vol;
		}
	}

	// body force
	for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
	{
		Particle_rb &pcl = md.pcls[md.bfxs[bf_id].pcl_id];
		double bf_tmp = pcl.vol * pcl.density * md.bfxs[bf_id].bf;
		pcl.pn1->fx_ext += pcl.N1 * bf_tmp;
		pcl.pn2->fx_ext += pcl.N2 * bf_tmp;
		pcl.pn3->fx_ext += pcl.N3 * bf_tmp;
		pcl.pn4->fx_ext += pcl.N4 * bf_tmp;
	}
	for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
	{
		Particle_rb &pcl = md.pcls[md.bfys[bf_id].pcl_id];
		double bf_tmp = pcl.vol * pcl.density * md.bfys[bf_id].bf;
		pcl.pn1->fy_ext += pcl.N1 * bf_tmp;
		pcl.pn2->fy_ext += pcl.N2 * bf_tmp;
		pcl.pn3->fy_ext += pcl.N3 * bf_tmp;
		pcl.pn4->fy_ext += pcl.N4 * bf_tmp;
	}

	// surface force
	for (size_t tf_id = 0; tf_id < md.tx_num; ++tf_id)
	{
		Particle_rb &pcl = md.pcls[md.txs[tf_id].pcl_id];
		double tf_tmp = md.txs[tf_id].t;
		pcl.pn1->fx_ext += pcl.N1 * tf_tmp;
		pcl.pn2->fx_ext += pcl.N2 * tf_tmp;
		pcl.pn3->fx_ext += pcl.N3 * tf_tmp;
		pcl.pn4->fx_ext += pcl.N4 * tf_tmp;
	}
	for (size_t tf_id = 0; tf_id < md.ty_num; ++tf_id)
	{
		Particle_rb &pcl = md.pcls[md.tys[tf_id].pcl_id];
		double tf_tmp = md.tys[tf_id].t;
		pcl.pn1->fy_ext += pcl.N1 * tf_tmp;
		pcl.pn2->fy_ext += pcl.N2 * tf_tmp;
		pcl.pn3->fy_ext += pcl.N3 * tf_tmp;
		pcl.pn4->fy_ext += pcl.N4 * tf_tmp;
	}

	// update nodal acceleration
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_rb &n = md.nodes[n_id];
		if (n.m != 0.0)
		{
			n.ax = (n.fx_ext - n.fx_int) / n.m;
			n.ay = (n.fy_ext - n.fy_int) / n.m;
		}
	}
	// acceleration boundary conditions
	for (size_t ax_id = 0; ax_id < md.ax_num; ++ax_id)
		md.nodes[md.axs[ax_id].node_id].ax = md.axs[ax_id].a;
	for (size_t ay_id = 0; ay_id < md.ay_num; ++ay_id)
		md.nodes[md.ays[ay_id].node_id].ax = md.ays[ay_id].a;

	self.dtime = self.max_dt;
	// adjust time step
	// particle shouldn't travel more than 1/10 element each steps
	double alpha_h = md.h * self.h_elem_raio;
	double v_max = 0.0, v_tmp;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_rb &n = md.nodes[n_id];
		v_tmp = (abs(n.vx) + sqrt(n.vx*n.vx + 4.0*abs(n.ax)*alpha_h)) / 2.0;
		if (v_max < v_tmp) v_max = v_tmp;
		v_tmp = (abs(n.vy) + sqrt(n.vy*n.vy + 4.0*abs(n.ay)*alpha_h)) / 2.0;
		if (v_max < v_tmp) v_max = v_tmp;
	}
	for (size_t vx_id = 0; vx_id < md.vx_num; ++vx_id)
	{
		v_tmp = md.vxs[vx_id].v;
		if (v_max < v_tmp) v_max = v_tmp;
	}
	for (size_t vy_id = 0; vy_id < md.vy_num; ++vy_id)
	{
		v_tmp = md.vys[vy_id].v;
		if (v_max < v_tmp) v_max = v_tmp;
	}
	if (v_max != 0.0)
	{
		double dt_tmp = alpha_h / v_max;
		//min_dt < dt_tmp < max_dt
		if (self.dtime > dt_tmp)
			self.dtime = dt_tmp;
	}

	// update nodal velocity
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_rb &n = md.nodes[n_id];
		if (n.m != 0.0)
		{
			n.vx /= n.m;
			n.vx += n.ax * self.dtime;
			n.vy /= n.m;
			n.vy += n.ay * self.dtime;
		}
	}
	// velocity boundary conditions
	for (size_t vx_id = 0; vx_id < md.vx_num; ++vx_id)
	{
		Node_rb &n = md.nodes[md.vxs[vx_id].node_id];
		n.vx = md.vxs[vx_id].v;
		n.ax = 0.0;
	}
	for (size_t vy_id = 0; vy_id < md.vy_num; ++vy_id)
	{
		Node_rb &n = md.nodes[md.vys[vy_id].node_id];
		n.vy = md.vys[vy_id].v;
		n.ay = 0.0;
	}

	md.rigid_body.predict_motion_from_ext_force(self.dtime);

	//// Contact detection
	double x1, y1, x2, y2, x3, y3, x4, y4;
	md.rigid_body.init_transformation();
	md.rigid_body.get_bounding_box(x1, y1, x2, y2, x3, y3, x4, y4, md.h);
	md.rasterize_rect_on_grid(x1, y1, x2, y2, x3, y3, x4, y4);
	double dist, nx, ny;
	double f_con, fx_con, fy_con;
	double K = 100.0; // Panelty factor, stiffness of interpenetration
	int res;
	for (size_t elem_id = 0; elem_id < md.elem_num; ++elem_id)
	{
		Element_rb &elem = md.elems[elem_id];
		if (elem.pcls && elem.has_rigid_object)
		{
			for (Particle_rb *pcl_iter= elem.pcls; pcl_iter; pcl_iter = pcl_iter->next_by_elem)
			{
				Particle_rb &pcl = *pcl_iter;
				Point pcl_pos = { pcl.x, pcl.y };
				res = md.rigid_body.distance_from_boundary(pcl_pos, dist, nx, ny, md.h);
				if (res < 0) continue; // not in contact
				// modify distance to consider size of material point
				dist += sqrt(pcl.vol / 4.0); // treated as square
				if (dist > 0) // is in contact
				{
					// normal force
					f_con = K * dist;
					fx_con = f_con * nx;
					fy_con = f_con * ny;
					// add force on rigid body
					md.rigid_body.add_con_force(-fx_con, -fy_con, pcl.x, pcl.y);
					// add force on material point to bg grid
					pcl.pn1->fx_con += fx_con * pcl.N1;
					pcl.pn1->fy_con += fy_con * pcl.N1;
					pcl.pn2->fx_con += fx_con * pcl.N2;
					pcl.pn2->fy_con += fy_con * pcl.N2;
					pcl.pn3->fx_con += fx_con * pcl.N3;
					pcl.pn3->fy_con += fy_con * pcl.N3;
					pcl.pn4->fx_con += fx_con * pcl.N4;
					pcl.pn4->fy_con += fy_con * pcl.N4;
				}
			}
		}
	}

	// corrected velocity and distance of nodes and rigid body
	double dt_square = self.dtime * self.dtime;
	double ax_corr, ay_corr;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_rb &n = md.nodes[n_id];
		if (n.fx_con != 0.0 && n.fy_con != 0.0)
		{
			ax_corr = n.fx_con / n.m;
			n.ax += ax_corr;
			n.vx += ax_corr * self.dtime;
		}
		if (n.fy_con != 0.0)
		{
			ay_corr = n.fy_con / n.m;
			n.ay += ay_corr;
			n.vy += ay_corr * self.dtime;
		}
	}
	// reapply velocity boundary conditions
	for (size_t vx_id = 0; vx_id < md.vx_num; ++vx_id)
	{
		Node_rb &n = md.nodes[md.vxs[vx_id].node_id];
		n.vx = md.vxs[vx_id].v;
		n.ax = 0.0;
	}
	for (size_t vy_id = 0; vy_id < md.vy_num; ++vy_id)
	{
		Node_rb &n = md.nodes[md.vys[vy_id].node_id];
		n.vy = md.vys[vy_id].v;
		n.ay = 0.0;
	}

	md.rigid_body.correct_motion_by_con_force(self.dtime);

	// map variables back to particles
	double de11, de12, de22, dw12, de_vol;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle_rb &pcl = md.pcls[pcl_id];
		if (pcl.pelem)
		{
			Node_rb &n1 = *pcl.pn1;
			Node_rb &n2 = *pcl.pn2;
			Node_rb &n3 = *pcl.pn3;
			Node_rb &n4 = *pcl.pn4;

			pcl.vx += (pcl.N1 * n1.ax + pcl.N2 * n2.ax
					 + pcl.N3 * n3.ax + pcl.N4 * n4.ax) * self.dtime;
			pcl.vy += (pcl.N1 * n1.ay + pcl.N2 * n2.ay
					 + pcl.N3 * n3.ay + pcl.N4 * n4.ay) * self.dtime;

			pcl.ux += (pcl.N1 * n1.vx + pcl.N2 * n2.vx
					 + pcl.N3 * n3.vx + pcl.N4 * n4.vx) * self.dtime;
			pcl.uy += (pcl.N1 * n1.vy + pcl.N2 * n2.vy
					 + pcl.N3 * n3.vy + pcl.N4 * n4.vy) * self.dtime;
			pcl.x = pcl.x_ori + pcl.ux;
			pcl.y = pcl.y_ori + pcl.uy;
			//std::cout << pcl.x_ori << ", " << pcl.y_ori << "\n";

			de11 = (pcl.dN1_dx * n1.vx + pcl.dN2_dx * n2.vx
				  + pcl.dN3_dx * n3.vx + pcl.dN4_dx * n4.vx) * self.dtime;
			de22 = (pcl.dN1_dy * n1.vy + pcl.dN2_dy * n2.vy
				  + pcl.dN3_dy * n3.vy + pcl.dN4_dy * n4.vy) * self.dtime;
			de12 = (pcl.dN1_dy * n1.vx + pcl.dN2_dy * n2.vx
				  + pcl.dN3_dy * n3.vx + pcl.dN4_dy * n4.vx
				  + pcl.dN1_dx * n1.vy + pcl.dN2_dx * n2.vy
				  + pcl.dN3_dx * n3.vy + pcl.dN4_dx * n4.vy) * 0.5 * self.dtime;
			dw12 = (pcl.dN1_dy * n1.vx + pcl.dN2_dy * n2.vx
				  + pcl.dN3_dy * n3.vx + pcl.dN4_dy * n4.vx
				  - pcl.dN1_dx * n1.vy - pcl.dN2_dx * n2.vy
				  - pcl.dN3_dx * n3.vy - pcl.dN4_dx * n4.vy) * 0.5 * self.dtime;

			// update strain (also assume that strain increment is Jaumann rate)
			//ppcl->de11 +=  ppcl->dw12 * ppcl->e12 * 2.0;
			//ppcl->de22 += -ppcl->dw12 * ppcl->e12 * 2.0;
			//ppcl->de12 +=  ppcl->dw12 * (ppcl->e22 - ppcl->e11);
			pcl.e11 += de11;
			pcl.e22 += de22;
			pcl.e12 += de12;

			// update stress
			double E_tmp = pcl.E / (1.0 + pcl.niu) / (1.0 - 2.0 * pcl.niu);
			double ds11, ds12, ds22;
			ds11 = E_tmp * ((1.0 - pcl.niu) * de11 + pcl.niu * de22);
			ds12 = 2.0 * pcl.E / (2.0 * (1.0 + pcl.niu)) * de12;
			ds22 = E_tmp * (pcl.niu * de11 + (1.0 - pcl.niu) * de22);

			/* ------------------------------------------------------------------
			Rotate as Jaumann rate:
			tensor_rate = tensor_Jaumann_rate + tensor * dW_T + dW * tensor
			------------------------------------------------------------------- */
			//ds11 +=  ppcl->dw12 * ppcl->s12 * 2.0;
			//ds22 += -ppcl->dw12 * ppcl->s12 * 2.0;
			//ds12 +=  ppcl->dw12 * (ppcl->s22 - ppcl->s11);
			pcl.s11 += ds11;
			pcl.s22 += ds22;
			pcl.s12 += ds12;

			de_vol = de11 + de22;
			pcl.density /= (1.0 + de_vol);
		}
	}

	return 0;
}
