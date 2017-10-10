#include "SpringForce.h"
#include <GL/glut.h>

#include <algorithm>	// to use find(InputIterator first, InputIterator last, const T& val)
#include <iterator>
#include <vector>

GravityForce::GravityForce( Vec3f gravity_constant ) :
  gravity( gravity_constant ) 
  {
  	is_spring = false;
  }
  
Vec3f GravityForce::force()
{
  return gravity;	// for unit mass
}

SpringForce::SpringForce(Particle *p1, Particle * p2, double dist, double ks, double kd) :
  m_p1(p1), m_p2(p2), m_dist(dist), m_ks(ks), m_kd(kd) {
  	is_spring = true;
  }
  
void SpringForce::update_index( std::vector<Particle*> pVector )
{
  for( index_p1 = 0; pVector[index_p1] != m_p1; index_p1++ )
  	;
  for( index_p1 = 0; pVector[index_p2] != m_p2; index_p2++ )
  	; 
} 

int SpringForce::index_of_p1()
{
  return index_p1; 
}

int SpringForce::index_of_p2()
{
  return index_p2;
}

const float air_res = 5.0f;

Vec3f SpringForce::force_on_p1()
{
  Vec3f delta_Position = m_p1->m_Position - m_p2->m_Position;
  Vec3f delta_Velocity = m_p1->m_Velocity - m_p2->m_Velocity;
  float delta_Distance = norm(delta_Position);
  
  Vec3f force_on_p1 = - ( delta_Position / delta_Distance ) * ( m_ks * ( delta_Distance - m_dist ) + m_kd * ( Dot( delta_Velocity, delta_Position / delta_Distance ) ) )
  					   - Dot( m_p1->m_Velocity, m_p1->m_Velocity ) * air_res * ( m_p1->m_Velocity );
  
  return force_on_p1;
}

Vec3f SpringForce::force_on_p2()
{

  Vec3f delta_Position = m_p2->m_Position - m_p1->m_Position;
  Vec3f delta_Velocity = m_p2->m_Velocity - m_p1->m_Velocity;
  float delta_Distance = norm(delta_Position);
  
  Vec3f force_on_p2 = - ( delta_Position / delta_Distance ) * ( m_ks * ( delta_Distance - m_dist ) + m_kd * ( Dot( delta_Velocity, delta_Position / delta_Distance ) ) )
  					  - ( Dot( m_p2->m_Velocity, m_p1->m_Velocity ) ) * air_res * ( m_p2->m_Velocity ); 
  
  return force_on_p2;
  
  // return ( -force_on_p1() );
}

void SpringForce::draw()
{
  glBegin( GL_LINES );
  glColor3f(0.6, 0.7, 0.8);
  glVertex2f( m_p1->m_Position[0], m_p1->m_Position[1] );
  glColor3f(0.6, 0.7, 0.8);
  glVertex2f( m_p2->m_Position[0], m_p2->m_Position[1] );
  glEnd();
}
