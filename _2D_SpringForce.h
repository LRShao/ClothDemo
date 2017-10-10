#pragma once

#include "Particle.h"
#include <vector>

inline float Dot( const Vec2f &u, const Vec2f &v )
{
  return (u[0]*v[0]) + (u[1]*v[1]);
}

class NonconstraintForce {
	public:
		bool is_spring;
		
		virtual void draw() {
			return;// just an interface for derived class SpringForce, does nothing
		}
};

class GravityForce: public NonconstraintForce {
	public:
		GravityForce();

		Vec2f force();
	
	private:
		Vec2f const gravity; 	
};

class SpringForce: public NonconstraintForce {
	public:
  		SpringForce(Particle *p1, Particle * p2, double dist, double ks, double kd);

  		void draw();
  		
  		void update_index( std::vector<Particle*> pVector );
  		int index_of_p1();
  		int index_of_p2();
  		
  		Vec2f force_on_p1();
  		Vec2f force_on_p2();

 	private:
  		Particle * const m_p1;   // particle 1
  		Particle * const m_p2;   // particle 2 
  		double const m_dist;     // rest length
  		double const m_ks, m_kd; // spring strength constants ( first as stiffness( coefficient in Hooke's law ), second as damping coeffecient )
  		int index_p1, index_p2;
};
