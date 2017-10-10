#include "Particle.h"
#include "SpringForce.h"	// definition of NonconstraintForce

#include <vector>
#include <cstdio>
#include <iostream>

#define DAMP 0.98f
#define RAND (((rand()%2000)/1000.f)-1.f)

std::vector<Vec2f> variableVector;		// in our case, should be [x1,v1,x2,v2,...,xn,vn]
std::vector<Vec2f> derivativeVector;	// should be [v1,f1/m1,...,vn,fn/mn]

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
	Vec2f zero_force(0.0, 0.0);
	
	for(ii=0; ii<size; ii++)
	{
		derivativeVector.push_back( pVector[ii]->m_Velocity );
		derivativeVector.push_back( zero_force );	// initialize force accumulator as zero
	}
	
	// Nonconstraint forces
	int fi, forceVectorSize = pNonconstraintForceVector.size();	// forceVectorSize = 2 * size
	for( fi=0; fi<forceVectorSize; fi++ )
	{
		printf ("\npNonconstrantForceVectorLoop, fi = %d\n", fi);
		if( pNonconstraintForceVector[fi]->is_spring )
		{
			printf ( "\n\n Running the Spring branch in nonconstraintforce:\n\n" );
			SpringForce* pCurrentForce = ( SpringForce* )( pNonconstraintForceVector[fi] );
			
			pCurrentForce->update_index( pVector );
			
			int index_of_p1 = pCurrentForce->index_of_p1(),
				index_of_p2 = pCurrentForce->index_of_p2();
				
			derivativeVector[ 2 * index_of_p1 + 1 ] += pCurrentForce->force_on_p1();	// Suppose that every particle has unit mass
			derivativeVector[ 2 * index_of_p2 + 1 ] += pCurrentForce->force_on_p2();
		}
		else
		{
			printf ( "\n\n Running the Gravity branch in nonconstraintforce:\n\n" );
			// either spring force or gravity force, should be applied to every particle
			for(ii=0; ii<size; ii++)
			{
				//GravityForce* pCurrentForce = ( GravityForce* )( pNonconstraintForceVector[fi] );
				// TEST
				Vec2f gravity(0,-0.02);
				derivativeVector[ 2 * ii + 1 ] += gravity;
				//derivativeVector[ 2 * ii + 1 ] += pCurrentForce->force();
			}						
		}
	}
	
	// TODO just for test, let us keep particle0 fixed, i.e. set its acceleration to zero
	derivativeVector[1] = zero_force;
	
	/* Unnecessary there, using default mass as 1g per particle */
	// Calculate acceleration from force and mass
	/*
	for(ii=0; ii<size; ii++)
	{
		if( pNonconstraintForceVector[ 2 * ii + 1 ] != zero_force );
			pNonconstraintForceVector[ 2 * ii + 1 ] /= pVector[ii]->mass;
	}
	*/	
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

void simulation_step( std::vector<Particle*> pVector, std::vector<NonconstraintForce*> pNonconstraintForceVector, float dt )
{
	/*****Euler Method******/
	/***********************
	// printf("\n==============================\nSIMULATION_STEP\n\nI.Getting variable data\n");	
	interface_get_variable_data( pVector );	
	
	int vi, size = variableVector.size();
	
	/*
	printf("variable vector:" );
	for(vi=0; vi<size; vi++)
	{
		std::cout << variableVector[vi] << "\t";
	}
	printf("\n");
	*
	
	interface_derivative_evaluation( pVector, pNonconstraintForceVector );
	/*
	printf("derivative vector:" );
	for(vi=0; vi<size; vi++)
	{
		std::cout << derivativeVector[vi] << "\t";
	}
	printf("\n");
	
	// int vi, size = variableVector.size();
	
	printf("\nII.Euler Step\n");
	*
	for(vi=0; vi<size; vi++)
	{
		variableVector[vi] += dt * derivativeVector[vi];
	}
	/*
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
	*
	
	interface_return_variable_data( pVector );
	
	variableVector.clear();
	derivativeVector.clear();
	*******************/
	
	/*****The Midpoint Method or Runge-Kutta 2 Method******
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
	
	*************************/
	
	/*****The Midpoint Method or Runge-Kutta 4 Method******/	
	interface_get_variable_data( pVector );
	std::vector<Vec2f> tempVariableVector = variableVector,
					   tempSumVector = variableVector;	
	
	int vi, size = variableVector.size();
	
	std::vector<Particle*> pTempParticleVector;
	
	for (int ii=0; ii<pVector.size(); ii++ )
	{
		pTempParticleVector.push_back( new Particle( Vec2f(0.0, 0.0) ) );	// Generate an independent copy of pVector called pTempParticleVector
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
	
	interface_return_variable_data( pVector );
	
	variableVector.clear();
	derivativeVector.clear();
}

