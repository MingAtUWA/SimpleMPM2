#ifndef __Disp_Con_Rigid_Circle_H__
#define __Disp_Con_Rigid_Circle_H__

class TriangleMeshToParticles;
class Model_T2D_CHM_s;

// Rigid ciricle to simulate T-bar penetration
// Displacement controlled
// Take counter-clockwise as positive

class DispConRigidCircle
{
	friend Model_T2D_CHM_s;

public:
	struct Particle
	{
		// position relative to the centroid
		double xr, yr;
		double vol; // volume
		// shape function
		double N1, N2, N3;
		double x, y;
		double vx, vy;
	};

	DispConRigidCircle() : pcl_num(0), pcls(nullptr) {}
	~DispConRigidCircle() { clear_pcl(); }
	
	void clear_pcl(void)
	{
		if (pcls)
		{
			delete[] pcls;
			pcls = nullptr;
		}
		pcl_num = 0;
	}

	int init(double _r, double _x, double _y, double max_pcl_size);
	void del_pcls_in_circle(TriangleMeshToParticles &tri2pcl, double exp_r = 0.0);
	void set_velocity(double _vx, double _vy, double _w);
	
protected:
	inline bool is_in_circle(double x, double y) noexcept
	{
		double x_diff = x - cen_x;
		double y_diff = y - cen_y;
		if (x_diff * x_diff + y_diff * y_diff > r2)
			return false;
		return true;
	}

	// update rigid body state
	void update(double dt);

	inline void add_reaction_force(double x, double y, double fx, double fy)
	{
		rfx += fx;
		rfy += fy;
		rm += (x - cen_x) * fy + (y - cen_y) * fx;
	}

public: // rigid body state
	struct State
	{
		double r, r2;
		double cen_x, cen_y, theta;
		double vx, vy, w;
		double rfx, rfy, rm;
		size_t pcl_num;
		Particle *pcls;
	};

protected:
	union
	{
		struct
		{
			// radius
			double r, r2;
			// position
			double cen_x, cen_y, theta;
			// velocity
			double vx, vy, w;
			// reaction force
			double rfx, rfy, rm;
			size_t pcl_num;
			Particle *pcls;
		};
		State state;
	};

public: // get result
	inline const State &get_state(void) noexcept { return state; }
	inline size_t get_pcl_num(void) noexcept { return pcl_num; }
	inline const Particle *get_pcls(void) noexcept { return pcls; }
	inline bool is_init(void) noexcept { return pcl_num != 0; }
};

#endif