#include "SimulationCore_pcp.h"

#include "Step_S2D_ME_s_RigidBody.h"
#include "TimeHistory_S2D_ME_s_RigidBody.h"

#include "Step_S2D_ME_s_RigidBody_Fric.h"
#include "TimeHistory_S2D_ME_s_RigidBody_Fric.h"

#include "ResultFile_XML.h"

ResultFile_XML::ResultFile_XML() {}

ResultFile_XML::~ResultFile_XML() { finalize(); }

int ResultFile_XML::init(const char *file_name)
{
	file.open(file_name, std::ios::binary | std::ios::out);
	if (!file.is_open())
		return -1;

	// xml version info
	const char *ver_info = "<?xml version=\"1.0\" encoding=\"ascii\"?>\n";
	file.write(ver_info, strlen(ver_info));

	return 0;
}

void ResultFile_XML::finalize(void)
{
	file.close();
}

// MPM - Rigid body contact problem
// Output model data
int ResultFile_XML::output(Model_S2D_ME_s_RigidBody &model)
{
	char str_buffer[512];

	// mesh
	const char *mesh_info = ""
		"<BackGroundMesh type = \"S2D\">\n"
		"    <h> %16.10e </h>\n"
		"    <x0> %16.10e </x0>\n"
		"    <xn> %16.10e </xn>\n"
		"    <y0> %16.10e </y0>\n"
		"    <yn> %16.10e </yn>\n"
		"    <elem_x_num> %zu </elem_x_num>\n"
		"    <elem_y_num> %zu </elem_y_num>\n"
		"</BackGroundMesh>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mesh_info, model.h,
			 model.x0, model.xn, model.y0, model.yn,
			 model.elem_x_num, model.elem_y_num);
	file.write(str_buffer, strlen(str_buffer));

	// rigid body
	TriangleMesh &rb_mesh = model.rigid_body.mesh;
	size_t node_num = rb_mesh.get_node_num();
	size_t elem_num = rb_mesh.get_elem_num();
	const char *rigid_obj_info = ""
		"<RigidObject type = \"TriangleMesh\">\n"
		"    <node_num> %zu </node_num>\n"
		"    <elem_num> %zu </elem_num>\n"
		"    <x_mc> %16.10e </x_mc>\n"
		"    <y_mc> %16.10e </y_mc>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), rigid_obj_info,
		node_num, elem_num, rb_mesh.get_x_mc(), rb_mesh.get_y_mc());
	file.write(str_buffer, strlen(str_buffer));
	// node coordinates
	file << "    <nodes>\n"
			"    <!-- index, x, y -->\n";
	const TriangleMesh::Node *rb_mesh_nodes = rb_mesh.get_nodes();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]),
				 "        %zu, %16.10e, %16.10e\n", n_id, 
				 rb_mesh_nodes[n_id].x, rb_mesh_nodes[n_id].y);
		file.write(str_buffer, strlen(str_buffer));
	}
	file << "    </nodes>\n";
	// element topology
	file << "    <elements>\n"
			"    <!-- index, node1, node2, node3 -->\n";
	const TriangleMesh::Element *rb_mesh_elems = rb_mesh.get_elems();
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		file << "        " << e_id
			 << ", " << rb_mesh_elems[e_id].n1
			 << ", " << rb_mesh_elems[e_id].n2
			 << ", " << rb_mesh_elems[e_id].n3 << "\n";
	}
	file << "    </elements>\n";
	// ending
	file << "</RigidObject>\n";

	// material point object
	const char *mp_obj_info = ""
		"<MaterialPointObject type = \"ME_2D\">\n"
		"    <pcl_num> %zu </pcl_num>\n"
		"</MaterialPointObject>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mp_obj_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));

	return 0;
}

// Output time history
int rf_out_func_imp_th_MPM_RigidBody_XML(TimeHistory &_th, ResultFile &_rf)
{
	TimeHistory_S2D_ME_s_RigidBody &th = static_cast<TimeHistory_S2D_ME_s_RigidBody &>(_th);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(_rf);
	std::fstream &file = rf.file;
	char str_buffer[512];

	// time history
	Step_S2D_ME_s_RigidBody &step = static_cast<Step_S2D_ME_s_RigidBody &>(th.get_step());
	const char *time_history_info = ""
		"<TimeHistory>\n"
		"    <substep_num> %zu </substep_num>\n"
		"    <total_substep_num> %zu </total_substep_num>\n"
		"    <current_time> %16.10e </current_time>\n"
		"    <total_time> %16.10e </total_time>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), time_history_info,
			 step.get_substep_num(), step.get_total_substep_num(),
			 step.get_current_time(), step.get_total_time());
	file.write(str_buffer, strlen(str_buffer));

	// rigid object
	Model_S2D_ME_s_RigidBody &model = static_cast<Model_S2D_ME_s_RigidBody &>(th.get_model());
	RigidBody &rb = model.rigid_body;
	const char *rigid_object_info = ""
		"    <RigidObject>\n"
		"        <x> %16.10e </x>\n"
		"        <y> %16.10e </y>\n"
		"        <theta> %16.10e </theta>\n"
		"        <vx> %16.10e </current_time>\n"
		"        <vy> %16.10e </total_time>\n"
		"        <vtheta> %16.10e </vtheta>\n"
		"        <Fx_contact> %16.10e </Fx_contact>\n"
		"        <Fy_contact> %16.10e </Fy_contact>\n"
		"        <M_contact> %16.10e </M_contact>\n"
		"    </RigidObject>\n";
	// resistence caused by contact
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), rigid_object_info,
			 rb.x, rb.y, rb.theta, rb.vx, rb.vy, rb.vtheta, rb.Fx_con, rb.Fy_con, rb.M_con);
	file.write(str_buffer, strlen(str_buffer));
	
	// output material points data
	const char *material_point_info = ""
		"    <MaterialPointObject>\n"
		"        <pcl_num> %zu </pcl_num>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), material_point_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));
	// field data: x, y, vol
	file << "        <field_data>\n"
		    "        <!-- x, y, vol -->\n";
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		file << "            " << model.pcls[pcl_id].x << ", " << model.pcls[pcl_id].y << ", "
			 << model.pcls[pcl_id].m / model.pcls[pcl_id].density << "\n";
	}
	file << "        </field_data>\n";
	// ending
	file << "    </MaterialPointObject>\n";

	// ending
	file << "</TimeHistory>\n";

	return 0;
}

const OutputFunc ResultFile_XML::out_func_th_MPM_RigidBody_XML = &rf_out_func_imp_th_MPM_RigidBody_XML;

int ResultFile_XML::output(Model_S2D_ME_s_RigidBody_Fric &model)
{
	char str_buffer[512];

	// mesh
	const char *mesh_info = ""
		"<BackGroundMesh type = \"S2D\">\n"
		"    <h> %16.10e </h>\n"
		"    <x0> %16.10e </x0>\n"
		"    <xn> %16.10e </xn>\n"
		"    <y0> %16.10e </y0>\n"
		"    <yn> %16.10e </yn>\n"
		"    <elem_x_num> %zu </elem_x_num>\n"
		"    <elem_y_num> %zu </elem_y_num>\n"
		"</BackGroundMesh>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mesh_info, model.h,
		model.x0, model.xn, model.y0, model.yn,
		model.elem_x_num, model.elem_y_num);
	file.write(str_buffer, strlen(str_buffer));

	// rigid body
	TriangleMesh &rb_mesh = model.rigid_body.mesh;
	size_t node_num = rb_mesh.get_node_num();
	size_t elem_num = rb_mesh.get_elem_num();
	const char *rigid_obj_info = ""
		"<RigidObject type = \"TriangleMesh\">\n"
		"    <node_num> %zu </node_num>\n"
		"    <elem_num> %zu </elem_num>\n"
		"    <x_mc> %16.10e </x_mc>\n"
		"    <y_mc> %16.10e </y_mc>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), rigid_obj_info,
		node_num, elem_num, rb_mesh.get_x_mc(), rb_mesh.get_y_mc());
	file.write(str_buffer, strlen(str_buffer));
	// node coordinates
	file << "    <nodes>\n"
		"    <!-- index, x, y -->\n";
	const TriangleMesh::Node *rb_mesh_nodes = rb_mesh.get_nodes();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]),
			"        %zu, %16.10e, %16.10e\n", n_id,
			rb_mesh_nodes[n_id].x, rb_mesh_nodes[n_id].y);
		file.write(str_buffer, strlen(str_buffer));
	}
	file << "    </nodes>\n";
	// element topology
	file << "    <elements>\n"
		"    <!-- index, node1, node2, node3 -->\n";
	const TriangleMesh::Element *rb_mesh_elems = rb_mesh.get_elems();
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		file << "        " << e_id
			<< ", " << rb_mesh_elems[e_id].n1
			<< ", " << rb_mesh_elems[e_id].n2
			<< ", " << rb_mesh_elems[e_id].n3 << "\n";
	}
	file << "    </elements>\n";
	// ending
	file << "</RigidObject>\n";

	// material point object
	const char *mp_obj_info = ""
		"<MaterialPointObject type = \"ME_2D\">\n"
		"    <pcl_num> %zu </pcl_num>\n"
		"</MaterialPointObject>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mp_obj_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));

	return 0;
}

// Output time history
int rf_out_func_imp_th_MPM_RigidBody_Fric_XML(TimeHistory &_th, ResultFile &_rf)
{
	TimeHistory_S2D_ME_s_RigidBody_Fric &th = static_cast<TimeHistory_S2D_ME_s_RigidBody_Fric &>(_th);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(_rf);
	std::fstream &file = rf.file;
	char str_buffer[512];

	// time history
	Step_S2D_ME_s_RigidBody_Fric &step = static_cast<Step_S2D_ME_s_RigidBody_Fric &>(th.get_step());
	const char *time_history_info = ""
		"<TimeHistory>\n"
		"    <substep_num> %zu </substep_num>\n"
		"    <total_substep_num> %zu </total_substep_num>\n"
		"    <current_time> %16.10e </current_time>\n"
		"    <total_time> %16.10e </total_time>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), time_history_info,
		step.get_substep_num(), step.get_total_substep_num(),
		step.get_current_time(), step.get_total_time());
	file.write(str_buffer, strlen(str_buffer));

	// rigid object
	Model_S2D_ME_s_RigidBody_Fric &model = static_cast<Model_S2D_ME_s_RigidBody_Fric &>(th.get_model());
	RigidBody &rb = model.rigid_body;
	const char *rigid_object_info = ""
		"    <RigidObject>\n"
		"        <x> %16.10e </x>\n"
		"        <y> %16.10e </y>\n"
		"        <theta> %16.10e </theta>\n"
		"        <vx> %16.10e </current_time>\n"
		"        <vy> %16.10e </total_time>\n"
		"        <vtheta> %16.10e </vtheta>\n"
		"        <Fx_contact> %16.10e </Fx_contact>\n"
		"        <Fy_contact> %16.10e </Fy_contact>\n"
		"        <M_contact> %16.10e </M_contact>\n"
		"    </RigidObject>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), rigid_object_info,
		rb.x, rb.y, rb.theta, rb.vx, rb.vy, rb.vtheta, rb.Fx_con, rb.Fy_con, rb.M_con);
	file.write(str_buffer, strlen(str_buffer));

	// output material points data
	const char *material_point_info = ""
		"    <MaterialPointObject>\n"
		"        <pcl_num> %zu </pcl_num>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), material_point_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));
	// field data: x, y, vol
	file << "        <field_data>\n"
		"        <!-- x, y, vol -->\n";
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		file << "            " << model.pcls[pcl_id].x << ", " << model.pcls[pcl_id].y << ", "
			<< model.pcls[pcl_id].m / model.pcls[pcl_id].density << "\n";
	}
	file << "        </field_data>\n";
	// ending
	file << "    </MaterialPointObject>\n";

	// ending
	file << "</TimeHistory>\n";

	return 0;
}

const OutputFunc ResultFile_XML::out_func_th_MPM_RigidBody_Fric_XML = &rf_out_func_imp_th_MPM_RigidBody_Fric_XML;
