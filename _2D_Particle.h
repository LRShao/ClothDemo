#pragma once

#include <gfx/vec2.h>
#include <gfx/vec3.h>

class Particle
{
public:

	Particle(const Vec2f & ConstructPos);
	virtual ~Particle(void);

	void reset();	// return the particle back to construction position and set velocity back to zero
	void draw();	// draw the particle as a square ( with white color and side length h = 0.03 )
	
	void operator = ( const Particle* pParticle );

	Vec2f m_ConstructPos;	// starting position
	Vec2f m_Position;	// current position
	Vec2f m_Velocity;	// current velocity
};
