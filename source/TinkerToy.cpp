// TinkerToy.cpp : Defines the entry point for the console application.

// Physics
#include "Particle.h"
#include "SpringForce.h"

// Extensible parts for constrained dynamics, unnecessary for cloth simulation
#include "RodConstraint.h"
#include "CircularWireConstraint.h"

// Screenshot
#include "imageio.h"
// Render
#include "shader.h"

// Graphics libraries
#include <GL/glew.h>
#include <GL/glut.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cstring>
#include <math.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace glm;

/* macros */
/* integration mode switch: 
 * Euler for Euler's method, accurate to O(dt)
 * Midpoint for midpoint method, accurate to O(dt^2)	
 * RK4 for Runge-Kutta4 method, accurate to O(dt^4)	
 */
const std::string MODE = "RK4";

/* external definitions (from solver.cpp) */
extern void simulation_step( std::vector<Particle*> pVector, std::vector<NonconstraintForce*> pNonconstraintForceVector, float dt, std::string mode );

/* global variables */

static int N;
static float dt, d;		// dt is time step in solver, ?d is the step size in dumping? 
static int dsim;		// if dsim == 0, simulation is in or will be set to initial state, else it is running.		
static int dump_frames;		// if is true, then frame dumping function is active
static int frame_number;	// the sequence number of current frame	

// static Particle *pList;
static std::vector<Particle*> pVector;	// keyword static means this variable can't be accessed from any other translation unit

static int win_id;		// window id returned by glutCreateWindow()
static int win_x, win_y;	// size of window
static int mouse_down[3];
static int mouse_release[3];
static int mouse_shiftclick[3];
static int omx, omy, mx, my;
static int hmx, hmy;

static std::vector<NonconstraintForce*> pNonconstraintForceVector;
static RodConstraint * delete_this_dummy_rod = NULL;
static CircularWireConstraint * delete_this_dummy_wire = NULL;

Vec3f rotate( Vec3f vector)		// identical to multiply by a rotation matrix
{
  // return Vec3f( sqrt(2.0)/2 * ( vector[0] - vector[1] ), -1, sqrt(2.0)/2 * ( vector[0] + vector[1] ));
  return Vec3f( -vector[1], -vector[2], vector[0] );
}

Vec3f rotate2( Vec3f vector)		// identical to multiply by a rotation matrix
{
  return Vec3f( vector[2], vector[1], vector[0] );
}

/*
----------------------------------------------------------------------
free/clear/allocate simulation data
----------------------------------------------------------------------
*/

static void free_data ( void )
{
	pVector.clear();
	if (delete_this_dummy_rod) {
		delete delete_this_dummy_rod;
		delete_this_dummy_rod = NULL;
	}
	if ( pNonconstraintForceVector.size() > 0 ) {
		pNonconstraintForceVector.clear();
	}
	if (delete_this_dummy_wire) {
		delete delete_this_dummy_wire;
		delete_this_dummy_wire = NULL;
	}
}

static void clear_data ( void )
{
	int ii, size = pVector.size();

	for(ii=0; ii<size; ii++){
		pVector[ii]->reset();
	}
}

static void init_system(void)
{
	/* Simple description of spring-particle model of cloth:
	   N X N particle grid, 
	   three type of springs
	   		stretch spring: grid side, connecting two nearest particles
	   		shear spring: diagonal of a grid block
			bend spring: contains two grid side, connecting every other two particles
	*/

	// define particle grid
	const double grid_length = 0.05;
	const double diagonal_length = sqrt(2.0) * grid_length;
	const int N = 20;			
	
	// offset from the center of the screen 
	const Vec3f x_positive_offset(grid_length, 0.0, 0.0 ),
				y_positive_offset(0.0, grid_length, 0.0);
			
	// Create particles as an N x N grid 
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			pVector.push_back(new Particle( rotate( i * y_positive_offset + j * x_positive_offset ) + Vec3f(0.5, 0.5, 0.0)  ) );
			
	// 1. Create universal gravity force
	// this force is just a dummy one
	pNonconstraintForceVector.push_back(new GravityForce( Vec3f(0, 0, 0) ));
	
	const float ks_stretch = 30, ks_shear = 30, ks_bend = 50,
				kd_stretch = 15, kd_shear = 15, kd_bend = 15;
	
	// original value, 60 and 20
	// 2. Stretch springs 
	// (N-1) x N horizontal spring constraint
	for(int i=0; i<N; i++)
		for(int j=0; j<(N-1); j++)
			pNonconstraintForceVector.push_back(new SpringForce(pVector[i*N+j], pVector[i*N+j+1], grid_length, ks_stretch, kd_stretch));
			
	// N x (N-1) vertical spring constraint
	for(int j=0; j<N; j++)
		for(int i=0; i<(N-1); i++)
			pNonconstraintForceVector.push_back(new SpringForce(pVector[i*N+j], pVector[(i+1)*N+j], grid_length, ks_stretch, kd_stretch));

	// 3. Shear springs 		
	// In total 2 X (N-1) x (N-1) diagonal spring constraint
	// those connecting bottom-left and upper-right particles
	for(int i=0; i<(N-1); i++)
		for(int j=0; j<(N-1); j++)
			pNonconstraintForceVector.push_back(new SpringForce(pVector[i*N+j], pVector[(i+1)*N+(j+1)], diagonal_length, ks_shear, kd_shear));

	// those connecting bottom-right and upper-left particles
	for(int i=1; i<N; i++)
		for(int j=0; j<(N-1); j++)
			pNonconstraintForceVector.push_back(new SpringForce(pVector[i*N+j], pVector[(i-1)*N+(j+1)], diagonal_length, ks_shear, kd_shear));

	// 4. Bend springs ( also a way to prevent self penetration in nearby region )
	for(int i=0; i<N; i++)
		for(int j=0; j<(N-2); j++)
			pNonconstraintForceVector.push_back(new SpringForce(pVector[i*N+j], pVector[i*N+j+2], 2 * grid_length, ks_bend, kd_bend));

	for(int j=0; j<N; j++)
		for(int i=0; i<(N-2); i++)
			pNonconstraintForceVector.push_back(new SpringForce(pVector[i*N+j], pVector[(i+2)*N+j], 2 * grid_length, ks_bend, kd_bend));	
}

/*
----------------------------------------------------------------------
OpenGL specific drawing routines
----------------------------------------------------------------------
*/

static void pre_display ( void )
{
	// display range is full window
	glViewport ( 0, 0, win_x, win_y );
	// Applies subsequent matrix operations to the projection matrix stack
	glMatrixMode ( GL_PROJECTION );
	// replace the current matrix with the identity matrix
	glLoadIdentity ();
	// define a 2D orthographic projection matrix. Four parameters are left, right, bottom, top clipping planes respectively
	gluOrtho2D ( -1.0, 1.0, -1.0, 1.0 );
	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_GREATER);
	// Specifies the current clearing values for the active color buffers	 
	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	// Clear all of the bound color buffers
	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}

static void post_display ( void )
{
	// Write frames if necessary.
	if (dump_frames) {
		const int FRAME_INTERVAL = 4;
		if ((frame_number % FRAME_INTERVAL) == 0) {
			const unsigned int w = glutGet(GLUT_WINDOW_WIDTH);
			const unsigned int h = glutGet(GLUT_WINDOW_HEIGHT);
			unsigned char * buffer = (unsigned char *) malloc(w * h * 4 * sizeof(unsigned char));
			if (!buffer)
				exit(-1);
			// glRasterPos2i(0, 0);
			glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
			char filename[13];
			sprintf(filename, "img%.5i.png", frame_number / FRAME_INTERVAL);
			printf("Dumped %s.\n", filename);
			saveImageRGBA(filename, buffer, w, h);
		}
	}
	frame_number++;
	
	// std::cout << "rendered frame: " << frame_number;
	
	glutSwapBuffers ();
}

static void draw_particles ( void )
{
	int size = pVector.size();

	for(int ii=0; ii< size; ii++)
	{
		pVector[ii]->draw();
	}
}

static void draw_forces ( void )
{
	// change this to iteration over full set
	std::vector<NonconstraintForce*>::iterator it;
	for( it = pNonconstraintForceVector.begin(); it < pNonconstraintForceVector.end(); it++ )
		if ( (*it)->is_spring )
			(*it)->draw();
}

static void draw_constraints ( void )
{
	// change this to iteration over full set
	if (delete_this_dummy_rod)
		delete_this_dummy_rod->draw();
	if (delete_this_dummy_wire)
		delete_this_dummy_wire->draw();
}

/*
----------------------------------------------------------------------
relates mouse movements to tinker toy construction
----------------------------------------------------------------------
*/

static void get_from_UI ()
{
	int i, j;
	// int size, flag;
	int hi, hj;
	// float x, y;
	if ( !mouse_down[0] && !mouse_down[2] && !mouse_release[0] 
	&& !mouse_shiftclick[0] && !mouse_shiftclick[2] ) return;

	i = (int)((       mx /(float)win_x)*N);
	j = (int)(((win_y-my)/(float)win_y)*N);

	if ( i<1 || i>N || j<1 || j>N ) return;

	if ( mouse_down[0] ) {

	}

	if ( mouse_down[2] ) {
	}

	hi = (int)((       hmx /(float)win_x)*N);
	hj = (int)(((win_y-hmy)/(float)win_y)*N);

	if( mouse_release[0] ) {
	}

	omx = mx;
	omy = my;
}

static void remap_GUI()
{
	int ii, size = pVector.size();
	for(ii=0; ii<size; ii++)
	{
		pVector[ii]->m_Position[0] = pVector[ii]->m_ConstructPos[0];
		pVector[ii]->m_Position[1] = pVector[ii]->m_ConstructPos[1];
	}
}

/*
----------------------------------------------------------------------
GLUT callback routines
----------------------------------------------------------------------
*/

static void key_func ( unsigned char key, int x, int y )
{
	switch ( key )
	{
	case 'c':
	case 'C':
		clear_data ();
		break;

	case 'd':
	case 'D':
		dump_frames = !dump_frames;
		break;

	case 'q':
	case 'Q':
		free_data ();
		exit ( 0 );
		break;

	case ' ':
		dsim = !dsim;
		break;
	}
}

static void mouse_func ( int button, int state, int x, int y )
{
	omx = mx = x;
	omx = my = y;

	if(!mouse_down[0]){hmx=x; hmy=y;}
	if(mouse_down[button]) mouse_release[button] = state == GLUT_UP;
	if(mouse_down[button]) mouse_shiftclick[button] = glutGetModifiers()==GLUT_ACTIVE_SHIFT;
	mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func ( int x, int y )	
{
	mx = x;
	my = y;
}

static void reshape_func ( int width, int height )
{
	glutSetWindow ( win_id );
	glutReshapeWindow ( width, height );

	win_x = width;
	win_y = height;
}

static void idle_func ( void )
{
	if ( dsim ) simulation_step( pVector, pNonconstraintForceVector, dt, MODE );
	else        {simulation_step( pVector, pNonconstraintForceVector, 0.0f, MODE );
				 /* get_from_UI();remap_GUI(); what is the purpose of this line ? */ }

	glutSetWindow ( win_id );
	glutPostRedisplay ();
}

// define the front and back color of the cloth
	const Vec3f FRONT_COLOR(1.0f, 0.647f, 0.0f);	// orange
	const Vec3f BACK_COLOR(0.0f, 1.0f, 0.498f); 	// spring green

// using vector function from <gfx/Vec3.h> 
Vec3f compute_lambertian_color( Vec3f p1_position, Vec3f p2_position, Vec3f p3_position ){
	Vec3f side1 = p2_position - p1_position;
	Vec3f side2 = p3_position - p1_position;
	
	Vec3f surface_normal;
	surface_normal = cross(side1, side2) / ( norm( cross(side1, side2) ) );
	
	if ( surface_normal[2] >= 0 )
		return FRONT_COLOR * surface_normal[2];
	else
		return BACK_COLOR * abs( surface_normal[2] );
}


static void display_func ( void )
{
	pre_display ();
	
	// ****************************************************************************
	GLuint VertexArrayID;
	
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders( "VertexShader.vertexshader", "FragmentShader.fragmentshader" );
	glBindAttribLocation(programID, 0, "vertexPosition_modelspace");
	glBindAttribLocation(programID, 1, "vertexColor");
	
	// Get a handle for our "MVP" uniform
	// Returns the index of the uniform variable name associated with the shader program
	// MVP stands for matrix-view-projection matrix
	GLuint MatrixID = glGetUniformLocation(programID, "MVP");

	/*
	// Projection matrix : 45?Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
	// TODO win_x = win_y = 512
	glm::mat4 Projection = glm::perspective(45.0f, 1.0f, 0.1f, 100.0f);
	// Camera matrix
	glm::mat4 View       = glm::lookAt(
								glm::vec3(0,0,0), // Camera is at (4,3,-3), in World Space
								glm::vec3(0,0,0), // and looks at the origin
								glm::vec3(0,0,0)  // Head is up (set to 0,-1,0 to look upside-down)
						   );
	// Model matrix : an identity matrix (model will be at the origin)
	glm::mat4 Model      = glm::mat4(1.0f);
	// Our ModelViewProjection : multiplication of our 3 matrices
	glm::mat4 MVP        = Projection * View * Model; // Remember, matrix multiplication is the other way around
	*/
	glm::mat4 MVP = glm::ortho(-1.0f,1.0f,-1.0f,1.0f,-1.0f,1.0f);

	const int N = 20;
		
	// Our vertices. Three consecutive floats give a 3D vertex; Three consecutive vertices give a triangle.
	// the cloth has N*N faces with 2 triangles each, so this makes 2*(N-1)*(N-1) triangles, and 2*(N-1)*(N-1)*3 vertices, and each vertex has 3 coordinate
	GLfloat g_vertex_buffer_data[ 2 * (N-1) * (N-1) * 3 * 3 ];
	for ( int i = 0; i < (N-1); i++ ){ 
		for ( int j = 0; j < (N-1); j++ ){
			// lower-left triangle
		
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) ] = pVector[ i * N + j ]->m_Position[0];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 1 ] = pVector[ i * N + j ]->m_Position[1];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 2 ] = pVector[ i * N + j ]->m_Position[2];
			
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 3 ] = pVector[ i * N + j + 1 ]->m_Position[0];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 4 ] = pVector[ i * N + j + 1 ]->m_Position[1];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 5 ] = pVector[ i * N + j + 1 ]->m_Position[2];
			
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 6 ] = pVector[ (i+1) * N + j ]->m_Position[0];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 7 ] = pVector[ (i+1) * N + j ]->m_Position[1];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 8 ] = pVector[ (i+1) * N + j ]->m_Position[2];
			
						
			// upper-right triangle
				
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 9 ] = pVector[ i * N + j + 1 ]->m_Position[0];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 10 ] = pVector[ i * N + j + 1 ]->m_Position[1];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 11 ] = pVector[ i * N + j + 1 ]->m_Position[2];
			
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 12 ] = pVector[ (i+1) * N + j ]->m_Position[0];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 13 ] = pVector[ (i+1) * N + j ]->m_Position[1];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 14 ] = pVector[ (i+1) * N + j ]->m_Position[2];
			
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 15 ] = pVector[ (i+1) * N + (j+1) ]->m_Position[0];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 16 ] = pVector[ (i+1) * N + (j+1) ]->m_Position[1];
			g_vertex_buffer_data[ 18 * ( i * (N-1) + j ) + 17 ] = pVector[ (i+1) * N + (j+1) ]->m_Position[2];
		}
	}
	

	// One color for each vertex. They were generated by applying lambertian shading
	GLfloat g_color_buffer_data[ 2 * (N-1) * (N-1) * 3 * 3 + 9 ];
	for ( int i = 0; i < (N-1); i++ ){ 
		for ( int j = 0; j < (N-1); j++ ){
			// lower-left triangle
			Vec3f color1 = compute_lambertian_color( pVector[ i * N + j ]->m_Position, pVector[ i * N + j + 1 ]->m_Position, pVector[ (i+1) * N + j ]->m_Position );
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) ] = color1[0];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 1 ] = color1[1];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 2 ] = color1[2];
			
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 3 ] = color1[0];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 4 ] = color1[1];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 5 ] = color1[2];
			
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 6 ] = color1[0];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 7 ] = color1[1];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 8 ] = color1[2];
						
			// upper-right triangle
			Vec3f color2 = compute_lambertian_color( pVector[ (i+1) * N + j ]->m_Position, pVector[ i * N + j + 1 ]->m_Position, pVector[ (i+1) * N + (j+1) ]->m_Position );
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 9 ] = color2[0];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 10 ] = color2[1];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 11 ] = color2[2];
			
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 12 ] = color2[0];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 13 ] = color2[1];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 14 ] = color2[2];
			
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 15 ] = color2[0];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 16 ] = color2[1];
			g_color_buffer_data[ 18 * ( i * (N-1) + j ) + 17 ] = color2[2];
		}
	}

	GLuint vertexbuffer;
	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_DYNAMIC_DRAW);

	GLuint colorbuffer;
	glGenBuffers(1, &colorbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_color_buffer_data), g_color_buffer_data, GL_DYNAMIC_DRAW);

		// Use our shader
		glUseProgram(programID);

		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glVertexAttribPointer(
			0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// 2nd attribute buffer : colors
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
		glVertexAttribPointer(
			1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
			3,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		
		// Draw the triangle !
		glDrawArrays(GL_TRIANGLES, 0, 2*(N-1)*(N-1)*3 ); // indices starting at 0 -> 2(N-1)*(N-1) triangles
		
		// glDrawArrays(GL_TRIANGLES, 0, 12*3); // 12*3 indices starting at 0 -> 12 triangles

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);


	// Cleanup VBO and shader
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteBuffers(1, &colorbuffer);
	glDeleteProgram(programID);
	glDeleteVertexArrays(1, &VertexArrayID);

	
	/*********************************************************/

	// Circular wire constraint : green
	// Spring force constraint : grey blue
	// Rod constraint: chocolate

	/*
	for( int i = 0; i < pVector.size(); i++ ){			
				printf("\npVector[%d]->m_Position =(%f, %f, %f)\n",
					    i, pVector[i]->m_Position[0], pVector[i]->m_Position[1], pVector[i]->m_Position[2]);
	}
	*/
	/*
	draw_forces();
	draw_constraints();
	draw_particles();
	*/

	/*********************************************************/	

	post_display ();
}


/*
----------------------------------------------------------------------
open_glut_window --- open a glut compatible window and set callbacks
----------------------------------------------------------------------
*/

static void open_glut_window ( void )
{	
	glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

	glutInitWindowPosition ( 0, 0 );
	glutInitWindowSize ( win_x, win_y );
	win_id = glutCreateWindow ( "ClothDemo" );
	
	glewExperimental = GL_TRUE; 
	glewInit();

	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	
	clear_data ();

	pre_display ();

	glutKeyboardFunc ( key_func );
	glutMouseFunc ( mouse_func );
	glutMotionFunc ( motion_func );
	glutReshapeFunc ( reshape_func );
	glutIdleFunc ( idle_func );
	glutDisplayFunc ( display_func );
}


/*
----------------------------------------------------------------------
main --- main routine
----------------------------------------------------------------------
*/

int main ( int argc, char ** argv )
{
	glutInit ( &argc, argv );

	if ( argc == 1 ) {
		N = 64;
		/*dt = 0.15f;*/
		dt = 0.015f;
		d = 5.f;
		fprintf ( stderr, "Using defaults : N=%d dt=%g d=%g\n",
			N, dt, d );
	} else {
		N = atoi(argv[1]);
		dt = atof(argv[2]);
		d = atof(argv[3]);
	}

	printf ( "\n\nHow to use this application:\n\n" );
	printf ( "\t Toggle construction/simulation display with the spacebar key\n" );
	printf ( "\t Dump frames by pressing the 'd' key\n" );
	printf ( "\t Quit by pressing the 'q' key\n" );

	dsim = 0;
	dump_frames = 0;
	frame_number = 0;
	
	init_system();
	
	win_x = 512;
	win_y = 512;
	open_glut_window ();

	glutMainLoop ();

	exit ( 0 );
}

