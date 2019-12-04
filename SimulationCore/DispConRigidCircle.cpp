#include "SimulationCore_pcp.h"

#include <cmath>
#include "TriangleMeshToParticles.h"

#include "DispConRigidCircle.h"

#define TWO_PI (3.14159265359 * 2.0)

int DispConRigidCircle::init(double _r, double _x, double _y, double max_pcl_size)
{
	r = _r;
	r2 = r * r;
	cen_x = _x;
	cen_y = _y;
	theta = 0.0;
	rfx = 0.0;
	rfy = 0.0;
	rm  = 0.0;
	// generate particles
	clear_pcl();
	double pcl_angle, pcl_size, pcl_r;
	pcl_num = size_t(ceil(TWO_PI * r / max_pcl_size));
	pcls = new Particle[pcl_num];
	pcl_angle = TWO_PI / double(pcl_num);
	pcl_size = pcl_angle * r;
	pcl_r = r - 0.5 * pcl_size;
	double angle;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		angle = double(pcl_id) * pcl_angle;
		pcl.xr = pcl_r * cos(angle);
		pcl.yr = pcl_r * sin(angle);
		pcl.vol = 0.5 * pcl_angle * (r2 - (r - pcl_size) * (r - pcl_size));
		pcl.x = cen_x + pcl.xr * cos(theta) + pcl.yr * -sin(theta);
		pcl.y = cen_y + pcl.xr * sin(theta) + pcl.yr *  cos(theta);
	}
	return 0;
}

void DispConRigidCircle::set_velocity(double _vx, double _vy, double _w)
{
	vx = _vx; vy = _vy; w = _w;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.vx = vx + pcl.yr * -w;
		pcl.vy = vy + pcl.xr *  w;
	}
}

void DispConRigidCircle::del_pcls_in_circle(TriangleMeshToParticles &tri2pcl)
{
	TriangleMeshToParticles::Particle *ppcl;
	ppcl = tri2pcl.first();
	while (tri2pcl.not_end_yet(ppcl))
	{
		if (is_in_circle(ppcl->x, ppcl->y))
		{
			TriangleMeshToParticles::Particle *ppcl_tmp;
			ppcl_tmp = ppcl;
			ppcl = tri2pcl.next(ppcl);
			tri2pcl.del_pcl(*ppcl_tmp);
			continue;
		}
		ppcl = tri2pcl.next(ppcl);
	}
}

void DispConRigidCircle::update(double dt)
{
	double dx, dy, dtheta;
	dx = vx * dt;
	dy = vy * dt;
	dtheta = w * dt;
	// update position
	cen_x += dx;
	cen_y += dy;
	theta += dtheta;
	// update particles position and velocity
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.x = cen_x + pcl.xr * cos(theta) + pcl.yr * -sin(theta);
		pcl.y = cen_y + pcl.xr * sin(theta) + pcl.yr *  cos(theta);
		pcl.vx = vx + pcl.yr * -w;
		pcl.vy = vy + pcl.xr *  w;
	}
}
