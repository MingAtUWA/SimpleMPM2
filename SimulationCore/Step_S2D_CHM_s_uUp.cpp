#include "SimulationCore_pcp.h"

#include "Step_S2D_CHM_s_uUp.h"

Step_S2D_CHM_s_uUp::Step_S2D_CHM_s_uUp() :
	Step(&solve_substep_S2D_CHM_s_uUp), model(nullptr) {}

Step_S2D_CHM_s_uUp::~Step_S2D_CHM_s_uUp() {}

int Step_S2D_CHM_s_uUp::init_calculation(void)
{
	if (is_first_step)
	{
		for (size_t pcl_id = 0; pcl_id < model->pcl_num; ++pcl_id)
		{
			Particle &pcl = model->pcls[pcl_id];
			pcl.pe = model->elems;
		}
	}

	for (size_t pcl_id = 0; pcl_id < model->pcl_num; ++pcl_id)
	{
		Particle &pcl = model->pcls[pcl_id];
		pcl.ux_s = 0.0;
		pcl.uy_s = 0.0;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux_f = 0.0;
		pcl.uy_f = 0.0;
	}
	
	// for debug
	out_file.open("debug_mat_out_CHM_s_uUp.csv", std::ios::binary | std::ios::out);

	map_from_pcl_to_node();

	return 0;
}

int Step_S2D_CHM_s_uUp::finalize_calculation(void)
{
	// for debug
	out_file.close();

	return 0;
}

namespace
{
typedef Model_S2D_CHM_s_uUp::Particle Particle_mpm;
typedef Model_S2D_CHM_s_uUp::Element Element_mpm;
typedef Model_S2D_CHM_s_uUp::Node Node_mpm;
typedef Model_S2D_CHM_s_uUp::DOF DOF;
};

int solve_substep_S2D_CHM_s_uUp(void *_self)
{
	Step_S2D_CHM_s_uUp &self = *(Step_S2D_CHM_s_uUp *)(_self);
	
	self.time_marching();
	self.form_global_residual_force();

	//self.map_from_pcl_to_node();

	//self.requilibration();
	//self.form_global_residual_force();

	return 0;
}

void Step_S2D_CHM_s_uUp::map_from_pcl_to_node(void)
{
	Model_S2D_CHM_s_uUp &md = *model;

	// init nodes
	size_t node_num = md.node_num;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		n.g_id = node_num;
		n.vol = 0.0;
		// solid phase
		n.ax_s = 0.0;
		n.ay_s = 0.0;
		n.vx_s = 0.0;
		n.vy_s = 0.0;
		// fluid phase
		n.ax_f = 0.0;
		n.ay_f = 0.0;
		n.vx_f = 0.0;
		n.vy_f = 0.0;
		n.p = 0.0;
	}
	// init elements
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
		md.elems[e_id].pcls = nullptr;

	// init particles
	double vol_N;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			if (md.init_pcl_shape_func(pcl))
			{
				pcl.vol = pcl.m_s / (pcl.density_s * (1.0 - pcl.n));
				Element_mpm &e = *pcl.pe;
				e.add_pcl(pcl);
				Node &n1 = md.nodes[e.n1_id];
				vol_N = pcl.vol * pcl.N1;
				n1.vol += vol_N;
				n1.p += pcl.p * vol_N;
				n1.ax_s += pcl.ax_s * vol_N;
				n1.ay_s += pcl.ay_s * vol_N;
				n1.vx_s += pcl.vx_s * vol_N;
				n1.vy_s += pcl.vy_s * vol_N;
				n1.ax_f += pcl.ax_f * vol_N;
				n1.ay_f += pcl.ay_f * vol_N;
				n1.vx_f += pcl.vx_f * vol_N;
				n1.vy_f += pcl.vy_f * vol_N;
				Node &n2 = md.nodes[e.n2_id];
				vol_N = pcl.vol * pcl.N2;
				n2.vol += vol_N;
				n2.p += pcl.p * vol_N;
				n2.ax_s += pcl.ax_s * vol_N;
				n2.ay_s += pcl.ay_s * vol_N;
				n2.vx_s += pcl.vx_s * vol_N;
				n2.vy_s += pcl.vy_s * vol_N;
				n2.ax_f += pcl.ax_f * vol_N;
				n2.ay_f += pcl.ay_f * vol_N;
				n2.vx_f += pcl.vx_f * vol_N;
				n2.vy_f += pcl.vy_f * vol_N;
				Node &n3 = md.nodes[e.n3_id];
				vol_N = pcl.vol * pcl.N3;
				n3.vol += vol_N;
				n3.p += pcl.p * vol_N;
				n3.ax_s += pcl.ax_s * vol_N;
				n3.ay_s += pcl.ay_s * vol_N;
				n3.vx_s += pcl.vx_s * vol_N;
				n3.vy_s += pcl.vy_s * vol_N;
				n3.ax_f += pcl.ax_f * vol_N;
				n3.ay_f += pcl.ay_f * vol_N;
				n3.vx_f += pcl.vx_f * vol_N;
				n3.vy_f += pcl.vy_f * vol_N;
				Node &n4 = md.nodes[e.n4_id];
				vol_N = pcl.vol * pcl.N4;
				n4.vol += vol_N;
				n4.p += pcl.p * vol_N;
				n4.ax_s += pcl.ax_s * vol_N;
				n4.ay_s += pcl.ay_s * vol_N;
				n4.vx_s += pcl.vx_s * vol_N;
				n4.vy_s += pcl.vy_s * vol_N;
				n4.ax_f += pcl.ax_f * vol_N;
				n4.ay_f += pcl.ay_f * vol_N;
				n4.vx_f += pcl.vx_f * vol_N;
				n4.vy_f += pcl.vy_f * vol_N;
			}
		}
	}

	cal_node_num = 0;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			Node &n1 = md.nodes[e.n1_id];
			if (n1.g_id == node_num)
			{
				n1.g_id = cal_node_num;
				++cal_node_num;
			}
			Node &n2 = md.nodes[e.n2_id];
			if (n2.g_id == node_num)
			{
				n2.g_id = cal_node_num;
				++cal_node_num;
			}
			Node &n4 = md.nodes[e.n4_id];
			if (n4.g_id == node_num)
			{
				n4.g_id = cal_node_num;
				++cal_node_num;
			}
			Node &n3 = md.nodes[e.n3_id];
			if (n3.g_id == node_num)
			{
				n3.g_id = cal_node_num;
				++cal_node_num;
			}
		}
	}

	dof_num = cal_node_num * 5; // ux_s, uy_s, ux_f, uy_f, p

	size_t *g_id_map = node_g_id_map_mem.resize(cal_node_num);
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		// cal nodal variables
		if (n.vol != 0.0)
		{
			n.ax_s /= n.vol;
			n.ay_s /= n.vol;
			n.vx_s /= n.vol;
			n.vy_s /= n.vol;
			n.ax_f /= n.vol;
			n.ay_f /= n.vol;
			n.vx_f /= n.vol;
			n.vy_f /= n.vol;
			n.p /= n.vol;
		}
		// form g_id_map
		if (n.g_id != node_num)
			g_id_map[n.g_id] = n_id;
	}

	//for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	//{
	//	Node &n = md.nodes[n_id];
	//	std::cout << n.g_id << " " 
	//			  << n.ax_s << " " << n.ay_s << " "
	//			  << n.vx_s << " " << n.vy_s << " "
	//			  << n.ax_f << " " << n.ay_f << " "
	//			  << n.vx_f << " " << n.vy_f << " "
	//			  << n.p << "\n";
	//}

	// apply inititial conditions, only support zero init cond
	for (size_t bc_id = 0; bc_id < md.usx_num; ++bc_id)
	{
		DisplacementBC &dbc = md.usxs[bc_id];
		Node &n = md.nodes[dbc.node_id];
		if (n.g_id != md.node_num)
		{
			n.ax_s = 0.0; // should = dbc.a(cur_t)
			n.vx_s = 0.0; // should = dbc.v(cur_t)
		}
	}
	for (size_t bc_id = 0; bc_id < md.usy_num; ++bc_id)
	{
		DisplacementBC &dbc = md.usys[bc_id];
		Node &n = md.nodes[dbc.node_id];
		if (n.g_id != md.node_num)
		{
			n.ay_s = 0.0; // should = dbc.a(cur_t)
			n.vy_s = 0.0; // should = dbc.v(cur_t)
		}
	}
	for (size_t bc_id = 0; bc_id < md.ufx_num; ++bc_id)
	{
		DisplacementBC &dbc = md.ufxs[bc_id];
		Node &n = md.nodes[dbc.node_id];
		if (n.g_id != md.node_num)
		{
			n.ax_f = 0.0; // should = dbc.a(cur_t)
			n.vx_f = 0.0; // should = dbc.v(cur_t)
		}
	}
	for (size_t bc_id = 0; bc_id < md.ufy_num; ++bc_id)
	{
		DisplacementBC &dbc = md.ufys[bc_id];
		Node &n = md.nodes[dbc.node_id];
		if (n.g_id != md.node_num)
		{
			n.ay_f = 0.0; // should = dbc.a(cur_t)
			n.vy_f = 0.0; // should = dbc.v(cur_t)
		}
	}
	for (size_t bc_id = 0; bc_id < md.pbc_num; ++bc_id)
	{
		PressureBC &pbc = md.pbcs[bc_id];
		Node &n = md.nodes[pbc.node_id];
		if (n.g_id != md.node_num)
			n.p = 0.0; // should = pbc.p(cur_t)
	}
}