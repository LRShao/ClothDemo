#include "Particle.h"
#include <GL/glut.h>

Particle::Particle(const Vec3f & ConstructPos) :
	m_ConstructPos(ConstructPos), m_Position(Vec3f(0.0, 0.0, 0.0)), m_Velocity(Vec3f(0.0, 0.0, 0.0)) {
}

Particle::~Particle(void) {
}

void Particle::reset() {
	m_Position = m_ConstructPos;
	m_Velocity = Vec3f(0.0, 0.0, 0.0);
}

void Particle::draw() {
	const double h = 0.03;
	glColor3f(1.f, 1.f, 1.f); 
	glBegin(GL_QUADS);
	glVertex2f(m_Position[0]-h/2.0, m_Position[1]-h/2.0);
	glVertex2f(m_Position[0]+h/2.0, m_Position[1]-h/2.0);
	glVertex2f(m_Position[0]+h/2.0, m_Position[1]+h/2.0);
	glVertex2f(m_Position[0]-h/2.0, m_Position[1]+h/2.0);
	glEnd();
}

void Particle::operator = ( const Particle* pParticle ) {
	m_ConstructPos = pParticle->m_ConstructPos;
	m_Position = pParticle->m_Position;
	m_Velocity = pParticle->m_Velocity;
}
