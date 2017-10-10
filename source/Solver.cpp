#include "Particle.h"
#include "SpringForce.h"	// definition of NonconstraintForce

#include <vector>
#include <cstdio>
#include <iostream>

#include <cmath>
#include <cstring>

std::vector<Vec3f> variableVector;		// in our case, should be [x1,v1,x2,v2,...,xn,vn]
std::vector<Vec3f> derivativeVector;	// should be [v1,f1/m1,...,vn,fn/mn]

const Vec3f ZERO_FORCE(0.0, 0.0, 0.0);
const int N = 20;

bool are_nodes_adjacent( int i, int j ){
	return ( abs( i-j ) == N ) || ( ( abs( i-j ) == 1 ) && !( ( (i+j+1) % (2*N) ) == 0 ) )
	       || ( ( abs	( (i%N) - (j%N) ) == 1 ) && ( abs( i - j ) == N ) );
}

void interface_get_variable_data( std::vector<Particle*> pVector )
{
	int ii, size = pVector.size();
	for(ii=0; ii<size; ii++)
	{
		variableVector.push_back( pVector[ii]->m_Position );
		variableVector.push_back( pVector[ii]->m_Velocity );
	}
}

// accumulate all forces on a particle and evaluate acceleration
void interface_derivative_evaluation( std::vector<Particle*> pVector, std::vector<NonconstraintForce*> pNonconstraintForceVector )
{
	int ii, size = pVector.size();
	
	for(ii=0; ii<size; ii++)
	{
		derivativeVector.push_back( pVector[ii]->m_Velocity );
		derivativeVector.push_back( ZERO_FORCE );	// initialize force accumulator as zero
	}
	
	// Nonconstraint forces
	int fi, forceVectorSize = pNonconstraintForceVector.size();	// forceVectorSize = 2 * size
	for( fi=0; fi<forceVectorSize; fi++ )
	{
		if( pNonconstraintForceVector[fi]->is_spring )
		{
			// spring force, add to connected particles' force accumulator
			SpringForce* pCurrentForce = ( SpringForce* )( pNonconstraintForceVector[fi] );
			
			pCurrentForce->update_index( pVector );
			
			int index_of_p1 = pCurrentForce->index_of_p1(),
				index_of_p2 = pCurrentForce->index_of_p2();
			
			// Suppose that every particle has unit mass	
			derivativeVector[ 2 * index_of_p1 + 1 ] += pCurrentForce->force_on_p1();	
			derivativeVector[ 2 * index_of_p2 + 1 ] += pCurrentForce->force_on_p2();
		}
		else
		{
			// gravity force, should be applied to every particle
			for(ii=0; ii<size; ii++)
			{
				//GravityForce* pCurrentForce = ( GravityForce* )( pNonconstraintForceVector[fi] );
				//derivativeVector[ 2 * ii + 1 ] += pCurrentForce->force();
				
				// to increase efficiency
				Vec3f gravity(0.0,-0.03,0.0);
				derivativeVector[ 2 * ii + 1 ] += gravity;
			}						
		}
	}
	
	// TODO just for test, let us keep particle0 fixed, i.e. set its acceleration to zero
	for (int i=0; i<N; i++)
		derivativeVector[ 1 + 2 * i * N ] = ZERO_FORCE;
}

void interface_return_variable_data( std::vector<Particle*> pVector )
{
	int ii, size = pVector.size();
	for(ii=0; ii<size; ii++)
	{
		pVector[ii]->m_Position = variableVector[ 2 * ii ];
		pVector[ii]->m_Velocity = variableVector[ 2 * ii + 1];
	}
}

void euler_method( std::vector<Particle*> pVector, std::vector<NonconstraintForce*> pNonconstraintForceVector, float dt ) {
	/*****Euler's Method******/
	// printf("\n==============================\nSIMULATION_STEP\n\nI.Getting variable data\n");	
	interface_get_variable_data( pVector );	
	
	int vi, size = variableVector.size();
	
	printf("variable vector:" );
	for(vi=0; vi<size; vi++)
	{
		std::cout << variableVector[vi] << "\t";
	}
	printf("\n");
	
	interface_derivative_evaluation( pVector, pNonconstraintForceVector );

	printf("derivative vector:" );
	for(vi=0; vi<size; vi++)
	{
		std::cout << derivativeVector[vi] << "\t";
	}
	printf("\n");
	
	// int vi, size = variableVector.size();
	
	printf("\nII.Euler Step\n");

	for(vi=0; vi<size; vi++)
	{
		variableVector[vi] += dt * derivativeVector[vi];
	}

	printf("\nIII.Returning variable data\n");
	
	printf("variable vector:" );
	for(vi=0; vi<size; vi++)
	{
		std::cout << variableVector[vi] << "\t";
	}
	printf("\n");
	
	printf("derivative vector:" );
	for(vi=0; vi<size; vi++)
	{
		std::cout << derivativeVector[vi] << "\t";
	}
	printf("\n==============================\n");
	
	interface_return_variable_data( pVector );
	
	variableVector.clear();
	derivativeVector.clear();

}

void midpoint_method( std::vector<Particle*> pVector, std::vector<NonconstraintForce*> pNonconstraintForceVector, float dt ) {
	/*****The Midpoint Method or Runge-Kutta 2 Method******/
	
	interface_get_variable_data( pVector );	

	int vi, size = variableVector.size();

	interface_derivative_evaluation( pVector, pNonconstraintForceVector );

	for(vi=0; vi<size; vi++)
	{
		variableVector[vi] += dt * derivativeVector[vi] / 2;
	}

	interface_return_variable_data( pVector );

	for(vi=0; vi<size; vi++)
	{
		variableVector[vi] -= dt * derivativeVector[vi] / 2;
	}

	derivativeVector.clear();
	
	interface_derivative_evaluation( pVector, pNonconstraintForceVector );

	for(vi=0; vi<size; vi++)
	{
		variableVector[vi] += dt * derivativeVector[vi];
	}

	interface_return_variable_data( pVector );

	variableVector.clear();
	derivativeVector.clear();
}

void runge_kutta4_method( std::vector<Particle*> pVector, std::vector<NonconstraintForce*> pNonconstraintForceVector, float dt ) {
	/*****Runge-Kutta 4 Method******/
	interface_get_variable_data( pVector );
	std::vector<Vec3f> tempVariableVector = variableVector,
					   tempSumVector = variableVector;	
	
	int vi, size = variableVector.size();
	
	// Generate an independent copy of pVector called pTempParticleVector
	std::vector<Particle*> pTempParticleVector;
	
	// deep copy construction
	for (int ii=0; ii<pVector.size(); ii++ )
	{
		pTempParticleVector.push_back( new Particle( Vec3f(0.0, 0.0, 0.0) ) );	
		*( pTempParticleVector.back() ) = pVector[ii];
	}
	
	// k_1 = hf( x_0 , t_0 )
	interface_derivative_evaluation( pVector, pNonconstraintForceVector );	// f(x_0)
	
	for(vi=0; vi<size; vi++)
	{
		tempVariableVector[vi] = dt * derivativeVector[vi];		// k_1
		tempSumVector[vi] += dt * tempVariableVector[vi] / 6;	// x_0 + k_1 / 6
	}
	
	for(vi=0; vi<size; vi++)
	{
		variableVector[vi] += tempVariableVector[vi] / 2;
	}
	
	interface_return_variable_data( pVector );	// x_0 + k_1 / 2
	
	variableVector.clear();
	derivativeVector.clear();
	
	// k_2 = hf( x_0 + k_1 / 2 , t_0 + h / 2 )	
	interface_derivative_evaluation( pVector, pNonconstraintForceVector );	
	
	for(vi=0; vi<size; vi++)
	{
		tempVariableVector[vi] = dt * derivativeVector[vi];	// k_2
		tempSumVector[vi] += dt * tempVariableVector[vi] / 3;	// x_0 + k_1 / 6 + k_2 / 3 
	}
	
	interface_get_variable_data( pTempParticleVector );	
	
	for(vi=0; vi<size; vi++)
	{
		variableVector[vi] += tempVariableVector[vi] / 2;
	}
	
	interface_return_variable_data( pVector );	// x_0 + k_2 / 2
	
	variableVector.clear();
	derivativeVector.clear();
	
	// k_3 = hf( x_0 + k_2 / 2 , t_0 + h / 2 )
	interface_derivative_evaluation( pVector, pNonconstraintForceVector );	
	
	for(vi=0; vi<size; vi++)
	{
		tempVariableVector[vi] = dt * derivativeVector[vi];	// k_3
		tempSumVector[vi] += dt * tempVariableVector[vi] / 3;	// x_0 + k_1 / 6 + k_2 / 3 + k_3 / 3
	}
	
	interface_get_variable_data( pTempParticleVector );
	
	for(vi=0; vi<size; vi++)
	{
		variableVector[vi] += tempVariableVector[vi];
	}
	
	interface_return_variable_data( pVector );	// x_0 + k_3
	
	variableVector.clear();
	derivativeVector.clear();
	
	// k4 = hf( x_0 + k_3, t_0 + h )	
	interface_derivative_evaluation( pVector, pNonconstraintForceVector );	
	
	interface_get_variable_data( pTempParticleVector );	
	
	for(vi=0; vi<size; vi++)
	{
		variableVector[vi] = dt * derivativeVector[vi];		// k_4
		tempSumVector[vi] += dt * derivativeVector[vi] / 6;	// tempSum = x_0 + k_1 / 6 + k_2 / 3 + k_3 / 3 + k_4 /6
		variableVector[vi] = tempSumVector[vi];				// variable = tempSum
	}
	
	/*
	// a very simple self-collision detection mechanism described in https://graphics.stanford.edu/~mdfisher/cloth.html
	// approximate each node particle as a marble with certain radius
	const float NODE_RADIUS = 0.0225;
	for (int ii=0; ii < pVector.size(); ii++ ){
		for (int jj=(ii+1); jj < pVector.size(); jj++ )
		{
			if( ( norm( variableVector[ 2 * ii ] - variableVector[ 2 * jj ] ) < ( 2 * NODE_RADIUS ) )
				&& ! are_nodes_adjacent( ii, jj ) ) 
				{
					variableVector[ 2 * ii ] = pVector[ii]->m_Position;
					variableVector[ 2 * jj ] = pVector[jj]->m_Position;
					
					/*
					variableVector[ 2 * ii + 1 ] = -variableVector[ 2 * ii + 1] * 0.3;
					variableVector[ 2 * jj + 1 ] = -variableVector[ 2 * ii + 1] * 0.3;
					*
					
					variableVector[ 2 * ii + 1 ] = Vec3f(0.0f,0.0f,0.0f);
					variableVector[ 2 * jj + 1 ] = Vec3f(0.0f,0.0f,0.0f);
					
					/*printf("\nmodified nodes: %d and %d, original distance: %f\n", ii, jj, norm( pVector[ii]->m_Position - pVector[jj]->m_Position ) );*		
				
				}
		}		
	}
	
	for (int i=0; i<N; i++)
		variableVector[ 1 + 2 * i * N ] = Vec3f(0.0f,0.0f,0.0f);
	
	interface_return_variable_data( pVector );
	*/
	
	variableVector.clear();
	derivativeVector.clear();

}

void simulation_step( std::vector<Particle*> pVector, std::vector<NonconstraintForce*> pNonconstraintForceVector, float dt, std::string mode )
{
	if ( mode == "Euler" )
		euler_method( pVector, pNonconstraintForceVector, dt );
	else if ( mode == "Midpoint" )
		midpoint_method( pVector, pNonconstraintForceVector, dt );
	else if ( mode == "RK4")
		runge_kutta4_method( pVector, pNonconstraintForceVector, dt );
	else
		std::cout << "No matching integration mode!";	
}


