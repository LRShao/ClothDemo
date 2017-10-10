#pragma once

#include "Particle.h"
#include <vector>

inline float Dot( const Vec3f &u, const Vec3f &v )
{
  return (u[0]*v[0]) + (u[1]*v[1]) + (u[2]*v[2]);
}

class NonconstraintForce {
	public:
		bool is_spring;		// identify whether this force is gravity or spring force
		
		virtual void draw() {
			return;
		}					// a virtual function, if this is spring force, draw it, otherwise does nothing
};

class GravityForce: public NonconstraintForce {
	public:
		GravityForce( Vec3f gravity_constant );

		Vec3f force();
	
	private:
		Vec3f const gravity; 	
};

class SpringForce: public NonconstraintForce {
	public:
  		SpringForce(Particle *p1, Particle * p2, double dist, double ks, double kd);

  		void draw();
  		
  		void update_index( std::vector<Particle*> pVector );
  		int index_of_p1();
  		int index_of_p2();
  		
  		Vec3f force_on_p1();
  		Vec3f force_on_p2();

 	private:
  		Particle * const m_p1;   // particle 1
  		Particle * const m_p2;   // particle 2 
  		double const m_dist;     // rest length
  		double const m_ks, m_kd; // spring strength constants ( first as stiffness( coefficient in Hooke's law ), second as damping coeffecient )
  		int index_p1, index_p2;
};
