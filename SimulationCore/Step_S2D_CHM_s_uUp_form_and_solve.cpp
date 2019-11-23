#include "SimulationCore_pcp.h"

#define Keep_Newmark_Coefficients
#define Keep_Debug_Mat_And_Vec_Output
#include "Step_S2D_CHM_s_uUp.h"

void Step_S2D_CHM_s_uUp::form_and_solve_problem(Eigen::VectorXd &g_du_vec)
{
	Model_S2D_CHM_s_uUp &md = *model;
	Eigen::VectorXd g_fvec(dof_num);
	g_fvec.setZero();
	g_kmat_coefs.init(dof_num);
	
	size_t node_g_id, l2g_id_map[20];
	double e_kmat[20][20], e_fvec[20];	
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			// form elemental stiffness matrix and force vector
			form_elem_stiffness_mat_and_force_vec(e, e_kmat, e_fvec);

			// form map from global id to local id
			node_g_id = md.nodes[e.n1_id].g_id;
			l2g_id_map[0]  = n_id_to_dof_id(node_g_id, DOF::ux_s);
			l2g_id_map[4]  = n_id_to_dof_id(node_g_id, DOF::uy_s);
			l2g_id_map[8]  = n_id_to_dof_id(node_g_id, DOF::ux_f);
			l2g_id_map[12] = n_id_to_dof_id(node_g_id, DOF::uy_f);
			l2g_id_map[16] = n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = md.nodes[e.n2_id].g_id;
			l2g_id_map[1]  = n_id_to_dof_id(node_g_id, DOF::ux_s);
			l2g_id_map[5]  = n_id_to_dof_id(node_g_id, DOF::uy_s);
			l2g_id_map[9]  = n_id_to_dof_id(node_g_id, DOF::ux_f);
			l2g_id_map[13] = n_id_to_dof_id(node_g_id, DOF::uy_f);
			l2g_id_map[17] = n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = md.nodes[e.n3_id].g_id;
			l2g_id_map[2]  = n_id_to_dof_id(node_g_id, DOF::ux_s);
			l2g_id_map[6]  = n_id_to_dof_id(node_g_id, DOF::uy_s);
			l2g_id_map[10] = n_id_to_dof_id(node_g_id, DOF::ux_f);
			l2g_id_map[14] = n_id_to_dof_id(node_g_id, DOF::uy_f);
			l2g_id_map[18] = n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = md.nodes[e.n4_id].g_id;
			l2g_id_map[3]  = n_id_to_dof_id(node_g_id, DOF::ux_s);
			l2g_id_map[7]  = n_id_to_dof_id(node_g_id, DOF::uy_s);
			l2g_id_map[11] = n_id_to_dof_id(node_g_id, DOF::ux_f);
			l2g_id_map[15] = n_id_to_dof_id(node_g_id, DOF::uy_f);
			l2g_id_map[19] = n_id_to_dof_id(node_g_id, DOF::p);

			// add to global matrix and vector
			for (size_t l_id1 = 0; l_id1 < 20; ++l_id1)
			{
				size_t g_id1 = l2g_id_map[l_id1];
				for (size_t l_id2 = 0; l_id2 < 20; ++l_id2)
				{
					size_t g_id2 = l2g_id_map[l_id2];
					g_kmat_coefs.add_coefficient(g_id1, g_id2, e_kmat[l_id1][l_id2]);
				}
				g_fvec[g_id1] += e_fvec[l_id1];
			}
		}
	}

	// apply external force
	// traction
	size_t dof_g_id;
	for (size_t t_id = 0; t_id < md.tx_num; ++t_id)
	{
		TractionBC_MPM &tx = md.txs[t_id];
		Particle &pcl = md.pcls[tx.pcl_id];
		if (pcl.pe)
		{
			Element &elem = *pcl.pe;
			// node 1
			node_g_id = md.nodes[elem.n1_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_s);
			g_fvec[dof_g_id] += pcl.N1 * tx.t;
			// node 2
			node_g_id = md.nodes[elem.n2_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_s);
			g_fvec[dof_g_id] += pcl.N2 * tx.t;
			// node 3
			node_g_id = md.nodes[elem.n3_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_s);
			g_fvec[dof_g_id] += pcl.N3 * tx.t;
			// node 4
			node_g_id = md.nodes[elem.n4_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_s);
			g_fvec[dof_g_id] += pcl.N4 * tx.t;
		}
	}
	for (size_t t_id = 0; t_id < md.ty_num; ++t_id)
	{
		TractionBC_MPM &ty = md.tys[t_id];
		Particle &pcl = md.pcls[ty.pcl_id];
		if (pcl.pe)
		{
			Element &elem = *pcl.pe;
			// node 1
			node_g_id = md.nodes[elem.n1_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_s);
			g_fvec[dof_g_id] += pcl.N1 * ty.t;
			// node 2
			node_g_id = md.nodes[elem.n2_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_s);
			g_fvec[dof_g_id] += pcl.N2 * ty.t;
			// node 3
			node_g_id = md.nodes[elem.n3_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_s);
			g_fvec[dof_g_id] += pcl.N3 * ty.t;
			// node 4
			node_g_id = md.nodes[elem.n4_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_s);
			g_fvec[dof_g_id] += pcl.N4 * ty.t;
		}
	}
	// body force to be finished...

	size_t node_num = md.node_num;
	// Apply displacement boundary condition
	double dig_term;
	double *kmat_col = kmat_col_mem.resize(dof_num);
	for (size_t bc_id = 0; bc_id < md.usx_num; ++bc_id)
	{
		DisplacementBC &dbc = md.usxs[bc_id];
		node_g_id = md.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_s);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * dbc.u;
			g_fvec[dof_g_id] = dig_term * dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < md.usy_num; ++bc_id)
	{
		DisplacementBC &dbc = md.usys[bc_id];
		node_g_id = md.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_s);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * dbc.u;
			g_fvec[dof_g_id] = dig_term * dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < md.ufx_num; ++bc_id)
	{
		DisplacementBC &dbc = md.ufxs[bc_id];
		node_g_id = md.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_f);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * dbc.u;
			g_fvec[dof_g_id] = dig_term * dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < md.ufy_num; ++bc_id)
	{
		DisplacementBC &dbc = md.ufys[bc_id];
		node_g_id = md.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_f);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * dbc.u;
			g_fvec[dof_g_id] = dig_term * dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < md.pbc_num; ++bc_id)
	{
		PressureBC &pbc = md.pbcs[bc_id];
		node_g_id = md.nodes[pbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::p);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * pbc.p;
			g_fvec[dof_g_id] = dig_term * pbc.p;
		}
	}

	// solve
	Eigen::SparseMatrix<double> g_kmat(dof_num, dof_num);
	g_kmat.setFromTriplets(g_kmat_coefs.begin(), g_kmat_coefs.end());
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver(g_kmat);
	g_du_vec = solver.solve(g_fvec);
	
	static size_t call_id = 0;
	++call_id;
	//print_sparse_mat(g_kmat, out_file);
	//std::cout << call_id << "\n" << g_fvec << "\n\n" << g_du_vec << "\n";

	// reapply disp bc
	for (size_t bc_id = 0; bc_id < md.usx_num; ++bc_id)
	{
		DisplacementBC &dbc = md.usxs[bc_id];
		node_g_id = md.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_s);
			g_du_vec[dof_g_id] = dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < md.usy_num; ++bc_id)
	{
		DisplacementBC &dbc = md.usys[bc_id];
		node_g_id = md.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_s);
			g_du_vec[dof_g_id] = dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < md.ufx_num; ++bc_id)
	{
		DisplacementBC &dbc = md.ufxs[bc_id];
		node_g_id = md.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_f);
			g_du_vec[dof_g_id] = dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < md.ufy_num; ++bc_id)
	{
		DisplacementBC &dbc = md.ufys[bc_id];
		node_g_id = md.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_f);
			g_du_vec[dof_g_id] = dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < md.pbc_num; ++bc_id)
	{
		PressureBC &pbc = md.pbcs[bc_id];
		node_g_id = md.nodes[pbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::p);
			g_du_vec[dof_g_id] = pbc.p;
		}
	}
}

int Step_S2D_CHM_s_uUp::requilibration(void)
{
	Model_S2D_CHM_s_uUp &md = *model;

	Eigen::VectorXd g_du_vec(dof_num);
	form_and_solve_problem(g_du_vec);

	// update nodal variables
	double dt = dtime, dt2 = dt * dt;
	size_t node_g_id;
	double dvx_s, dvy_s, dvx_f, dvy_f;
	double dax_s, day_s, dax_f, day_f;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		node_g_id = n.g_id;
		if (node_g_id != md.node_num)
		{
			n.dux_s = g_du_vec[n_id_to_dof_id(node_g_id, DOF::ux_s)];
			n.duy_s = g_du_vec[n_id_to_dof_id(node_g_id, DOF::uy_s)];
			n.dux_f = g_du_vec[n_id_to_dof_id(node_g_id, DOF::ux_f)];
			n.duy_f = g_du_vec[n_id_to_dof_id(node_g_id, DOF::uy_f)];
			n.dp    = g_du_vec[n_id_to_dof_id(node_g_id, DOF::p)];
			dvx_s = gamma / (beta * dt) * n.dux_s - gamma / beta * n.vx_s + dt * (1.0 - gamma / (2.0 * beta)) * n.ax_s;
			dvy_s = gamma / (beta * dt) * n.duy_s - gamma / beta * n.vy_s + dt * (1.0 - gamma / (2.0 * beta)) * n.ay_s;
			dvx_f = gamma / (beta * dt) * n.dux_f - gamma / beta * n.vx_f + dt * (1.0 - gamma / (2.0 * beta)) * n.ax_f;
			dvy_f = gamma / (beta * dt) * n.duy_f - gamma / beta * n.vy_f + dt * (1.0 - gamma / (2.0 * beta)) * n.ay_f;
			dax_s = n.dux_s / (beta * dt2) - n.vx_s / (beta * dt) - n.ax_s / (2.0 * beta);
			day_s = n.duy_s / (beta * dt2) - n.vy_s / (beta * dt) - n.ay_s / (2.0 * beta);
			dax_f = n.dux_f / (beta * dt2) - n.vx_f / (beta * dt) - n.ax_f / (2.0 * beta);
			day_f = n.duy_f / (beta * dt2) - n.vy_f / (beta * dt) - n.ay_f / (2.0 * beta);
			n.vx_s += dvx_s;
			n.vy_s += dvy_s;
			n.vx_f += dvx_f;
			n.vy_f += dvy_f;
			n.ax_s += dax_s;
			n.ay_s += day_s;
			n.ax_f += dax_f;
			n.ay_f += day_f;
			n.p    += n.dp;
		}
	}

	// Update particle variables
	double de11, de22, de12;
	double ds11, ds22, ds12, dp;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			Node &n1 = md.nodes[e.n1_id];
			Node &n2 = md.nodes[e.n2_id];
			Node &n3 = md.nodes[e.n3_id];
			Node &n4 = md.nodes[e.n4_id];

			// update particle variables
			for (Particle *pcl_iter = e.pcls; pcl_iter; pcl_iter = pcl_iter->next)
			{
				Particle &pcl = *pcl_iter;
				
				// velocity
				pcl.vx_s = n1.vx_s * pcl.N1 + n2.vx_s * pcl.N2 + n3.vx_s * pcl.N3 + n4.vx_s * pcl.N4;
				pcl.vy_s = n1.vy_s * pcl.N1 + n2.vy_s * pcl.N2 + n3.vy_s * pcl.N3 + n4.vy_s * pcl.N4;
				pcl.vx_f = n1.vx_f * pcl.N1 + n2.vx_f * pcl.N2 + n3.vx_f * pcl.N3 + n4.vx_f * pcl.N4;
				pcl.vy_f = n1.vy_f * pcl.N1 + n2.vy_f * pcl.N2 + n3.vy_f * pcl.N3 + n4.vy_f * pcl.N4;

				// acceleration
				pcl.ax_s = n1.ax_s * pcl.N1 + n2.ax_s * pcl.N2 + n3.ax_s * pcl.N3 + n4.ax_s * pcl.N4;
				pcl.ay_s = n1.ay_s * pcl.N1 + n2.ay_s * pcl.N2 + n3.ay_s * pcl.N3 + n4.ay_s * pcl.N4;
				pcl.ax_f = n1.ax_f * pcl.N1 + n2.ax_f * pcl.N2 + n3.ax_f * pcl.N3 + n4.ax_f * pcl.N4;
				pcl.ay_f = n1.ay_f * pcl.N1 + n2.ay_f * pcl.N2 + n3.ay_f * pcl.N3 + n4.ay_f * pcl.N4;
				
				// strain increment
				de11 = n1.dux_s * pcl.dN1_dx + n2.dux_s * pcl.dN2_dx
					 + n3.dux_s * pcl.dN3_dx + n4.dux_s * pcl.dN4_dx;
				de22 = n1.duy_s * pcl.dN1_dy + n2.duy_s * pcl.dN2_dy
					 + n3.duy_s * pcl.dN3_dy + n4.duy_s * pcl.dN4_dy;
				de12 = (n1.dux_s * pcl.dN1_dy + n2.dux_s * pcl.dN2_dy
					  + n3.dux_s * pcl.dN3_dy + n4.dux_s * pcl.dN4_dy
					  + n1.duy_s * pcl.dN1_dx + n2.duy_s * pcl.dN2_dx
					  + n3.duy_s * pcl.dN3_dx + n4.duy_s * pcl.dN4_dx) * 0.5;
				pcl.e11 += de11;
				pcl.e22 += de22;
				pcl.e12 += de12;

				// update stress
				double E_tmp = pcl.E / (1.0 + pcl.niu) / (1.0 - 2.0 * pcl.niu);
				ds11 = E_tmp * ((1.0 - pcl.niu) * de11 + pcl.niu * de22);
				ds22 = E_tmp * (pcl.niu * de11 + (1.0 - pcl.niu) * de22);
				ds12 = pcl.E / (2.0 * (1.0 + pcl.niu)) * 2.0 * de12;
				pcl.s11 += ds11;
				pcl.s22 += ds22;
				pcl.s12 += ds12;

				// pore pressure
				dp = n1.dp * pcl.N1 + n2.dp * pcl.N2 + n3.dp * pcl.N3 + n4.dp * pcl.N4;
				pcl.p += dp;
				//pcl.p = n1.p * pcl.N1 + n2.p * pcl.N2 + n3.p * pcl.N3 + n4.p + pcl.N4;
			}
		}
	}
	
	return 0;
}

int Step_S2D_CHM_s_uUp::time_marching(void)
{
	Model_S2D_CHM_s_uUp &md = *model;

	Eigen::VectorXd g_du_vec(dof_num);
	form_and_solve_problem(g_du_vec);

	// Update nodal variables
	double dt = dtime, dt2 = dtime * dtime;
	size_t node_g_id;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		node_g_id = n.g_id;
		if (node_g_id != md.node_num)
		{
			n.dux_s = g_du_vec[n_id_to_dof_id(node_g_id, DOF::ux_s)];
			n.duy_s = g_du_vec[n_id_to_dof_id(node_g_id, DOF::uy_s)];
			n.dux_f = g_du_vec[n_id_to_dof_id(node_g_id, DOF::ux_f)];
			n.duy_f = g_du_vec[n_id_to_dof_id(node_g_id, DOF::uy_f)];
			n.dp = g_du_vec[n_id_to_dof_id(node_g_id, DOF::p)];
			n.dvx_s = gamma / (beta * dt) * n.dux_s - gamma / beta * n.vx_s + dt * (1.0 - gamma / (2.0 * beta)) * n.ax_s;
			n.dvy_s = gamma / (beta * dt) * n.duy_s - gamma / beta * n.vy_s + dt * (1.0 - gamma / (2.0 * beta)) * n.ay_s;
			n.dvx_f = gamma / (beta * dt) * n.dux_f - gamma / beta * n.vx_f + dt * (1.0 - gamma / (2.0 * beta)) * n.ax_f;
			n.dvy_f = gamma / (beta * dt) * n.duy_f - gamma / beta * n.vy_f + dt * (1.0 - gamma / (2.0 * beta)) * n.ay_f;
			n.dax_s = n.dux_s / (beta * dt2) - n.vx_s / (beta * dt) - n.ax_s / (2.0 * beta);
			n.day_s = n.duy_s / (beta * dt2) - n.vy_s / (beta * dt) - n.ay_s / (2.0 * beta);
			n.dax_f = n.dux_f / (beta * dt2) - n.vx_f / (beta * dt) - n.ax_f / (2.0 * beta);
			n.day_f = n.duy_f / (beta * dt2) - n.vy_f / (beta * dt) - n.ay_f / (2.0 * beta);
			n.vx_s += n.dvx_s;
			n.vy_s += n.dvy_s;
			n.vx_f += n.dvx_f;
			n.vy_f += n.dvy_f;
			n.ax_s += n.dax_s;
			n.ay_s += n.day_s;
			n.ax_f += n.dax_f;
			n.ay_f += n.day_f;
			n.p += n.dp;

		}
	}

	// Update particle position and size
	double de11, de22, de12, de_vol_s;
	double ds11, ds22, ds12, dp;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			Node &n1 = md.nodes[e.n1_id];
			Node &n2 = md.nodes[e.n2_id];
			Node &n3 = md.nodes[e.n3_id];
			Node &n4 = md.nodes[e.n4_id];

			// update particle variables
			for (Particle *pcl_iter = e.pcls; pcl_iter; pcl_iter = pcl_iter->next)
			{
				Particle &pcl = *pcl_iter;

				// velocity
				pcl.vx_s += n1.dvx_s * pcl.N1 + n2.dvx_s * pcl.N2 + n3.dvx_s * pcl.N3 + n4.dvx_s * pcl.N4;
				pcl.vy_s += n1.dvy_s * pcl.N1 + n2.dvy_s * pcl.N2 + n3.dvy_s * pcl.N3 + n4.dvy_s * pcl.N4;
				pcl.vx_f += n1.dvx_f * pcl.N1 + n2.dvx_f * pcl.N2 + n3.dvx_f * pcl.N3 + n4.dvx_f * pcl.N4;
				pcl.vy_f += n1.dvy_f * pcl.N1 + n2.dvy_f * pcl.N2 + n3.dvy_f * pcl.N3 + n4.dvy_f * pcl.N4;

				// acceleration
				pcl.ax_s += n1.dax_s * pcl.N1 + n2.dax_s * pcl.N2 + n3.dax_s * pcl.N3 + n4.dax_s * pcl.N4;
				pcl.ay_s += n1.day_s * pcl.N1 + n2.day_s * pcl.N2 + n3.day_s * pcl.N3 + n4.day_s * pcl.N4;
				pcl.ax_f += n1.dax_f * pcl.N1 + n2.dax_f * pcl.N2 + n3.dax_f * pcl.N3 + n4.dax_f * pcl.N4;
				pcl.ay_f += n1.day_f * pcl.N1 + n2.day_f * pcl.N2 + n3.day_f * pcl.N3 + n4.day_f * pcl.N4;
				
				// displacement
				pcl.ux_s += n1.dux_s * pcl.N1 + n2.dux_s * pcl.N2 + n3.dux_s * pcl.N3 + n4.dux_s * pcl.N4;
				pcl.uy_s += n1.duy_s * pcl.N1 + n2.duy_s * pcl.N2 + n3.duy_s * pcl.N3 + n4.duy_s * pcl.N4;
				pcl.ux_f += n1.dux_f * pcl.N1 + n2.dux_f * pcl.N2 + n3.dux_f * pcl.N3 + n4.dux_f * pcl.N4;
				pcl.uy_f += n1.duy_f * pcl.N1 + n2.duy_f * pcl.N2 + n3.duy_f * pcl.N3 + n4.duy_f * pcl.N4;
				// position
				pcl.x = pcl.x_ori + pcl.ux_s;
				pcl.y = pcl.y_ori + pcl.uy_s;

				// strain increment
				de11 = n1.dux_s * pcl.dN1_dx + n2.dux_s * pcl.dN2_dx
					 + n3.dux_s * pcl.dN3_dx + n4.dux_s * pcl.dN4_dx;
				de22 = n1.duy_s * pcl.dN1_dy + n2.duy_s * pcl.dN2_dy
					 + n3.duy_s * pcl.dN3_dy + n4.duy_s * pcl.dN4_dy;
				de12 = (n1.dux_s * pcl.dN1_dy + n2.dux_s * pcl.dN2_dy
					  + n3.dux_s * pcl.dN3_dy + n4.dux_s * pcl.dN4_dy
					  + n1.duy_s * pcl.dN1_dx + n2.duy_s * pcl.dN2_dx
					  + n3.duy_s * pcl.dN3_dx + n4.duy_s * pcl.dN4_dx) * 0.5;
				pcl.e11 += de11;
				pcl.e22 += de22;
				pcl.e12 += de12;

				// update stress
				double E_tmp = pcl.E / (1.0 + pcl.niu) / (1.0 - 2.0 * pcl.niu);
				ds11 = E_tmp * ((1.0 - pcl.niu) * de11 + pcl.niu * de22);
				ds22 = E_tmp * (pcl.niu * de11 + (1.0 - pcl.niu) * de22);
				ds12 = pcl.E / (2.0 * (1.0 + pcl.niu)) * 2.0 * de12;
				pcl.s11 += ds11;
				pcl.s22 += ds22;
				pcl.s12 += ds12;

				// pore pressure
				dp = n1.dp * pcl.N1 + n2.dp * pcl.N2 + n3.dp * pcl.N3 + n4.dp * pcl.N4;
				//pcl.p += dp;
				pcl.p = n1.p * pcl.N1 + n2.p * pcl.N2 + n3.p * pcl.N3 + n4.p * pcl.N4;
				
				// porosity
				de_vol_s = de11 + de22;
				pcl.n = (de_vol_s + pcl.n) / (1.0 + de_vol_s);

				// density
				pcl.density_f /= (1.0 + dp / pcl.Kf);
			}
		}
	}

	//// Update nodal variables
	//size_t node_g_id;
	//for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	//{
	//	Node &n = md.nodes[n_id];
	//	node_g_id = n.g_id;
	//	if (node_g_id != md.node_num)
	//	{
	//		n.dux_s = g_du_vec[n_id_to_dof_id(node_g_id, DOF::ux_s)];
	//		n.duy_s = g_du_vec[n_id_to_dof_id(node_g_id, DOF::uy_s)];
	//		n.dux_f = g_du_vec[n_id_to_dof_id(node_g_id, DOF::ux_f)];
	//		n.duy_f = g_du_vec[n_id_to_dof_id(node_g_id, DOF::uy_f)];
	//		n.dp    = g_du_vec[n_id_to_dof_id(node_g_id, DOF::p)];
	//	}
	//}

	//// Update particle position and size
	//double de11, de22, de_vol_s, dp;
	//for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	//{
	//	Element &e = md.elems[e_id];
	//	if (e.pcls)
	//	{
	//		Node &n1 = md.nodes[e.n1_id];
	//		Node &n2 = md.nodes[e.n2_id];
	//		Node &n3 = md.nodes[e.n3_id];
	//		Node &n4 = md.nodes[e.n4_id];
	//		
	//		// update particle variables
	//		for (Particle *pcl_iter = e.pcls; pcl_iter; pcl_iter = pcl_iter->next)
	//		{
	//			Particle &pcl = *pcl_iter;

	//			// displacement
	//			pcl.ux_s += n1.dux_s * pcl.N1 + n2.dux_s * pcl.N2 + n3.dux_s * pcl.N3 + n4.dux_s * pcl.N4;
	//			pcl.uy_s += n1.duy_s * pcl.N1 + n2.duy_s * pcl.N2 + n3.duy_s * pcl.N3 + n4.duy_s * pcl.N4;
	//			pcl.ux_f += n1.dux_f * pcl.N1 + n2.dux_f * pcl.N2 + n3.dux_f * pcl.N3 + n4.dux_f * pcl.N4;
	//			pcl.uy_f += n1.duy_f * pcl.N1 + n2.duy_f * pcl.N2 + n3.duy_f * pcl.N3 + n4.duy_f * pcl.N4;
	//			// position
	//			pcl.x = pcl.x_ori + pcl.ux_s;
	//			pcl.y = pcl.y_ori + pcl.uy_s;

	//			// porosity
	//			de11 = n1.dux_s * pcl.dN1_dx + n2.dux_s * pcl.dN2_dx + n3.dux_s * pcl.dN3_dx + n4.dux_s * pcl.dN4_dx;
	//			de22 = n1.duy_s * pcl.dN1_dy + n2.duy_s * pcl.dN2_dy + n3.duy_s * pcl.dN3_dy + n4.duy_s * pcl.dN4_dy;
	//			de_vol_s = de11 + de22;
	//			pcl.n = (de_vol_s + pcl.n) / (1.0 + de_vol_s);

	//			// density
	//			dp = n1.dp * pcl.N1 + n2.dp * pcl.N2 + n3.dp * pcl.N3 + n4.dp * pcl.N4;
	//			pcl.density_f /= (1.0 + dp / pcl.Kf);
	//		}
	//	}
	//}
			
	return 0;
}

