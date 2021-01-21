//
//  mpm_integrator.hpp
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2020/12/18.
//  Copyright © 2020年 Kentaro Nagasawa. All rights reserved.
//

#ifndef mpm_integrator_hpp
#define mpm_integrator_hpp

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <Eigen/SVD>
#include "omp.h"

#include "mpm_state.h"
#ifdef MPM2D
#include "../2DInterface/mpm_setting.h"
#endif

#ifdef MPM3D
#include "../3DInterface/mpm_setting.h"
#include "../3DInterface/Situation_Setting.h"
extern MPM::SituationSets				situation_set;
extern MPM::BoundaryCondition	boundary_conditions;
#endif



using Vectord = Vector<double, SET::dim>;
using Vectori = Vector<int, SET::dim>;
using Matrixd = Matrix<double, SET::dim>;
using Matrixi = Matrix<int, SET::dim>;

using MatrixXd = MatrixX<double, SET::dim>;

namespace MPM {
	class Basisfunctions
	{
		//std::vector<Vectord> N;
		//Basisfunctions();
	public:
		std::vector<Vectord> N;
		void set(const Vectord& d_bg);
		void init();
		double tensor(const Vectori &node);
		VectorXd vec();
	};

	auto BasisFunction(const Vectord& d_bg);

	struct Integrater
	{
		MPM::Basisfunctions basis_f;
		std::vector<MPM::Basisfunctions> basis_fs;
		MPM::SimulationState sim_state;
		Integrater() = default;
		//Integrater(MPM::SimulationState& _sim_setting);
		//Integrater() = default;

		void set_integrater(MPM::SimulationState& _sim_setting);
		void integrate();
		
		void position_update();
		void boundary_colision();
		
		void particle2grid();
		void gridAdvance();
		void grid2particle();
		void advection();
		void deformation_update();
		
		void stateCheck();
		void updateState();
		
	};

	

}

inline std::vector<Vectord> Basisfunction(const Vectord &d_bg)
{
	std::vector<Vectord> N;
	N.push_back( 0.5*( 1.5*Vectord::Ones(SET::dim) - d_bg ).cwiseAbs2() );
	N.push_back( 0.75*Vectord::Ones(SET::dim) - ( d_bg - Vectord::Ones(SET::dim) ).cwiseAbs2() );
	N.push_back( 0.5*(d_bg - 0.5*Vectord::Ones(SET::dim) ).cwiseAbs2() );
	return N;
}


inline void MPM::Basisfunctions::set(const Vectord &d_bg)
{
	//N.push_back( 0.5*( 1.5*Vectord::Ones(SET::dim) - d_bg ).cwiseAbs2() );
	//N.push_back( 0.75*Vectord::Ones(SET::dim) - ( d_bg - Vectord::Ones(SET::dim) ).cwiseAbs2() );
	//N.push_back( 0.5*(d_bg - 0.5*Vectord::Ones(SET::dim) ).cwiseAbs2() );
	N[0] = ( 0.5*( 1.5*Vectord::Ones(SET::dim) - d_bg ).cwiseAbs2() );
	N[1] = ( 0.75*Vectord::Ones(SET::dim) - ( d_bg - Vectord::Ones(SET::dim) ).cwiseAbs2() );
	N[2] = ( 0.5*(d_bg - 0.5*Vectord::Ones(SET::dim) ).cwiseAbs2() );

}

inline void MPM::Basisfunctions::init()
{
	N.push_back( Vectord::Ones(SET::dim) );
	N.push_back( Vectord::Ones(SET::dim) );
	N.push_back( Vectord::Ones(SET::dim) );
	//std::cout << "N's size is " << N.size() << std::endl;
}

inline VectorXd MPM::Basisfunctions::vec()
{
	const int Max = std::pow(3,SET::dim);
	VectorXd vec = VectorXd::Zero(Max);
	const Vectori w = Vectori::Constant(3);
	for(int i=0;i<Max;++i)
	{
		vec(i) = tensor(Vectori{MPM::flat2node(i, w)});
	}
	return std::move(vec);
}




inline double MPM::Basisfunctions::tensor(const Vectori &node)
{
	double N_val = 1.0;
	for(int i=0; i<SET::dim; ++i)
	{
		N_val *= N[node(i)](i);
	}
	return  std::move(N_val);
}


inline void MPM::Integrater::set_integrater(MPM::SimulationState& _sim_setting)
{
	sim_state = _sim_setting;
	basis_f.init();
	
	#ifdef USEOPENMP
	//Set per_thread basis functions

	basis_fs.push_back(MPM::Basisfunctions{});

	for( int i = 0; i < omp_get_max_threads(); i++ )
	{
		basis_fs.push_back(MPM::Basisfunctions{});
	}
	
	for( int i = 0; i < omp_get_max_threads(); i++ )
	{
		basis_fs[i].init();
	}
	#endif
	
}

inline void MPM::Integrater::integrate()
{
	stateCheck();
	
	particle2grid();
	gridAdvance();
	grid2particle();
	deformation_update();
	advection();
	
	updateState();
}

inline void MPM::Integrater::position_update()
{
sim_state.mp.Position.row(1).array() -= 0.01;
}


void MPM::Integrater::particle2grid()
{
	const double width_inv = sim_state.SimulationGridsWidth_inv;
	
	//Initialize ValiedGrid
	for(const auto& i : sim_state.mg.ValiedGridNodeIndices){
			sim_state.mg.Velocity.col(i) = Vectord::Zero();
			sim_state.mg.Mass[i] = 0;
	}
	sim_state.mg.ValiedGridNodeIndices.clear();
	
	const int NP=sim_state.mp.numpoints;
	#ifndef USEOPENMP
	for(int i=0; i<NP; ++i)
	{
		const Matrixd FP = sim_state.c->cal_particle_stress(sim_state.mp.F[i], sim_state.mp.J[i],sim_state.mp.b[i]);
		
		const double p_vol = sim_state.mp.vol[i];
		const double particle_mass = sim_state.mp.m[i];
		
		const Matrixd Stress = - (sim_state.dt * p_vol) * (4 * width_inv * width_inv * FP);
		const Matrixd Affine = Stress + sim_state.alpha*particle_mass * sim_state.mp.C[i];
		//const Matrixd Affine = Stress;
		
		const Vectori base_grid =	((sim_state.mp.Position.col(i) - sim_state.SimulationBoxMin ) *width_inv - Vectord::Constant(0.5)).cast<int>();
		const Vectord d_from_basegrid = (sim_state.mp.Position.col(i) - sim_state.SimulationBoxMin )*width_inv - base_grid.cast<double>();
		
		basis_f.set(d_from_basegrid.cast<double>());

		const Vectord mv = sim_state.mp.Velocity.col(i) * particle_mass;
		const Vectori w = Vectori::Constant(3);

				for(int l=0, Max = std::pow(3,SET::dim); l<Max; ++l)
		{
			const Vectori node 			= MPM::flat2node(l, w);
			const Vectori base_node = base_grid + node;
			const Vectord distance_from_index = (node.cast<double>() - d_from_basegrid) * sim_state.SimulationGridsWidth;

			sim_state.mg.Velocity.col( MPM::node2flat(base_node, sim_state.SimulationGridsnum) ) += basis_f.tensor(node) * ( mv + Affine * distance_from_index );
			sim_state.mg.Mass[ MPM::node2flat(base_node, sim_state.SimulationGridsnum) ] += particle_mass * basis_f.tensor(node);
			
			sim_state.mg.ValiedGridNodeIndices.push_back(MPM::node2flat(base_node, sim_state.SimulationGridsnum));
		}
	}
	
	#else //IF use OpenMP
	
	//Set per_thread_memories
	static std::vector<MatrixXd> 							per_thread_Velocity;
	static std::vector<std::vector<double>> 	per_thread_mass;
	static std::vector<std::vector<int>>		 	per_thread_ValiedGridNodeIndices;
	
	if( per_thread_Velocity.empty() )
	{
		per_thread_Velocity.resize( omp_get_max_threads() );
	}
	if( per_thread_mass.empty() )
	{
		per_thread_mass.resize( omp_get_max_threads() );
	}

	if( per_thread_ValiedGridNodeIndices.empty() )
	{
		per_thread_ValiedGridNodeIndices.resize( omp_get_max_threads() );
	}

	assert( int(per_thread_Velocity.size()) == omp_get_max_threads() );
	assert( int(per_thread_mass.size()) == omp_get_max_threads() );
	assert( int(per_thread_ValiedGridNodeIndices.size()) == omp_get_max_threads() );
	
	#pragma omp parallel for
	for( std::vector<MatrixXd>::size_type idx = 0; idx < per_thread_Velocity.size(); idx++ )
	{
		per_thread_Velocity[idx].setZero( SET::dim,sim_state.mg.numgrids );
	}
	
	#pragma omp parallel for
	for( std::vector<std::vector<double>>::size_type idx = 0; idx < per_thread_mass.size(); idx++ )
	{
		per_thread_mass[idx].resize( sim_state.mg.numgrids );
		std::fill(per_thread_mass[idx].begin(), per_thread_mass[idx].end(), 0.0);
	}
	
	#pragma omp parallel for
	for( std::vector<std::vector<int>>::size_type idx = 0; idx < per_thread_ValiedGridNodeIndices.size(); idx++ )
	{
		per_thread_ValiedGridNodeIndices[idx].clear();
	}

	
	#pragma omp parallel for
	for(int i=0; i<NP; ++i)
	{
		//get omp thread ID
		const int tid{ omp_get_thread_num() };
		assert( tid >= 0 );

		const Matrixd FP = sim_state.c->cal_particle_stress(sim_state.mp.F[i], sim_state.mp.J[i],sim_state.mp.b[i]);
		
		const double p_vol = sim_state.mp.vol[i];
		const double particle_mass = sim_state.mp.m[i];
				
		const Matrixd Stress = - (sim_state.dt * p_vol) * (4 * width_inv * width_inv * FP);
		const Matrixd Affine = Stress + sim_state.alpha*particle_mass * sim_state.mp.C[i];
		//const Matrixd Affine = Stress;
		
		const Vectori base_grid =	((sim_state.mp.Position.col(i) - sim_state.SimulationBoxMin ) *width_inv - Vectord::Constant(0.5)).cast<int>();
		const Vectord d_from_basegrid = (sim_state.mp.Position.col(i) - sim_state.SimulationBoxMin )*width_inv - base_grid.cast<double>();
		
		//set basis function on omp thread ID
		basis_fs[tid].set(d_from_basegrid.cast<double>());


		const Vectord mv = sim_state.mp.Velocity.col(i) * particle_mass;
		const Vectori w = Vectori::Constant(3);

		for(int l=0, Max = std::pow(3,SET::dim); l<Max; ++l)
		{
			const Vectori node 			= MPM::flat2node(l, w);
			const Vectori base_node = base_grid + node;
			const Vectord distance_from_index = (node.cast<double>() - d_from_basegrid) * sim_state.SimulationGridsWidth;

			per_thread_Velocity[tid].col( MPM::node2flat(base_node, sim_state.SimulationGridsnum) ) += basis_fs[tid].tensor(node) * ( mv + Affine * distance_from_index );

			per_thread_mass[tid][ MPM::node2flat(base_node, sim_state.SimulationGridsnum) ] += particle_mass * basis_fs[tid].tensor(node);
			
			per_thread_ValiedGridNodeIndices[tid].push_back(MPM::node2flat(base_node, sim_state.SimulationGridsnum));
			//sim_state.mg.ValiedGridNodeIndices.push_back(MPM::node2flat(base_node, sim_state.SimulationGridsnum));

		}
	}
	
	//#pragma omp parallel for
	for( std::vector<MatrixXd>::size_type idx = 0; idx < per_thread_Velocity.size(); idx++ )
	{
		sim_state.mg.Velocity += per_thread_Velocity[idx];
	}
	
	#pragma omp parallel for
	for(int i=0; i< sim_state.mg.numgrids; ++i)
	{
		for( std::vector<std::vector<double>>::size_type idx = 0; idx < per_thread_mass.size(); idx++ )
		{
			sim_state.mg.Mass[i] += per_thread_mass[idx][i];
		}
	}

	for( std::vector<MatrixXd>::size_type idx = 0; idx < per_thread_ValiedGridNodeIndices.size(); idx++ )
	{
		sim_state.mg.ValiedGridNodeIndices.insert(sim_state.mg.ValiedGridNodeIndices.end(), per_thread_ValiedGridNodeIndices[idx].begin(), per_thread_ValiedGridNodeIndices[idx].end());
	}

	#endif

	
	
}

void MPM::Integrater::gridAdvance()
{
	
	//Remove duplicates in ValiedGridNode
	//Todo: Make this more quick
	std::sort(sim_state.mg.ValiedGridNodeIndices.begin(), sim_state.mg.ValiedGridNodeIndices.end());
	decltype(sim_state.mg.ValiedGridNodeIndices)::iterator result = std::unique(sim_state.mg.ValiedGridNodeIndices.begin(), sim_state.mg.ValiedGridNodeIndices.end());
	sim_state.mg.ValiedGridNodeIndices.erase(result, sim_state.mg.ValiedGridNodeIndices.end());
	

		//for(int i=0,NP=sim_state.mg.numgrids;i<NP;++i)
	const int Valied_num = sim_state.mg.ValiedGridNodeIndices.size();
	#ifdef USEOPENMP
	#pragma omp parallel for
	#endif
	for(int j=0; j < Valied_num; j++)
	{
			const int i = sim_state.mg.ValiedGridNodeIndices[j];
		//if(sim_state.mg.Mass[i] > 0)
		//{
			sim_state.mg.Velocity.col(i) = sim_state.mg.Velocity.col(i)/sim_state.mg.Mass[i];
			Vectord gravity = Vectord::Zero();
			gravity.y() = -981;
			sim_state.mg.Velocity.col(i) += sim_state.dt * gravity;
			
			// Boundary and Collision behaviours ------------------
			const double bd_bottom = sim_state.SimulationGridsWidth*0.0;
			if(sim_state.mg.Position.col(i).y() < ( bd_bottom + sim_state.SimulationBoxMinBound.y() ) )
			{
				sim_state.mg.Velocity.col(i) = Vectord::Constant(0.0);
				//sim_state.mg.Velocity.col(i).y() = std::max(0.0,sim_state.mg.Velocity.col(i).y());
			}

			if( ( sim_state.mg.Position.col(i).array() < ( sim_state.SimulationBoxMinBound +  bd_bottom*Vectord::Ones() ).array() ).any() || (sim_state.mg.Position.col(i).array() > (sim_state.SimulationBoxMaxBound - bd_bottom*Vectord::Ones() ).array() ).any() )
			{
				if(sim_state.mg.Position.col(i).y() < bd_bottom + sim_state.SimulationBoxMinBound.y() )
				{
					continue;
				}
				sim_state.mg.Velocity.col(i) = Vectord::Constant(0.0);
			}
			
			#ifdef MPM3D //Apply 3D rigid body collision
			sim_state.mg.Velocity.col(i) = situation_set.boundary_conditions.ComputeBoundaryVelocity(sim_state.mg.Position.col(i), sim_state.mg.Velocity.col(i));
			#endif//----------------------------------------------------
			
		//}
		
	}
}

void MPM::Integrater::grid2particle()
{
	const double width_inv = sim_state.SimulationGridsWidth_inv;

	sim_state.mp.Velocity = Eigen::MatrixXd::Zero(SET::dim,sim_state.mp.numpoints);
	std::fill(sim_state.mp.C.begin(), sim_state.mp.C.end(), Matrixd::Zero());

	
	

	const int NP=sim_state.mp.numpoints;
	#ifndef USEOPENMP
	for(int i=0; i<NP; ++i)
	{
		const Vectori base_grid =	((sim_state.mp.Position.col(i) - sim_state.SimulationBoxMin ) *width_inv - Vectord::Constant(0.5)).cast<int>();
		const Vectord d_from_basegrid = (sim_state.mp.Position.col(i) - sim_state.SimulationBoxMin )*width_inv - base_grid.cast<double>();
		
		
		basis_f.set(d_from_basegrid.cast<double>());
		
		const Vectori w = Vectori::Constant(3);
		const int Max = std::pow(3,SET::dim);
		for(int l=0; l<Max; ++l)
		{
			const Vectori node 			= MPM::flat2node(l, w);
			const Vectori base_node = base_grid + node;
			const Vectord distance_from_index = (node.cast<double>() - d_from_basegrid);

			sim_state.mp.Velocity.col(i)+= basis_f.tensor(node)*sim_state.mg.Velocity.col( MPM::node2flat(base_node, sim_state.SimulationGridsnum)  ) ;
			
			sim_state.mp.C[i] += 4 * width_inv * basis_f.tensor(node) * sim_state.mg.Velocity.col(MPM::node2flat(base_node, sim_state.SimulationGridsnum) )*(distance_from_index).transpose();
		}
	}
	#else //IF use OpenMP
	
	//Set per_thread memories
	static std::vector<MatrixXd> 							per_thread_Velocity;
	static std::vector<std::vector<Matrixd>> 	per_thread_C;
	if( per_thread_Velocity.empty() )
	{
		per_thread_Velocity.resize( omp_get_max_threads() );
	}

	if( per_thread_C.empty() )
	{
		per_thread_C.resize( omp_get_max_threads() );
	}
	
	assert( int(per_thread_Velocity.size()) == omp_get_max_threads() );
	assert( int(per_thread_C.size()) == omp_get_max_threads() );

	#pragma omp parallel for
	for( std::vector<MatrixXd>::size_type idx = 0; idx < per_thread_Velocity.size(); idx++ )
	{
		per_thread_Velocity[idx].setZero( SET::dim,sim_state.mp.numpoints );
	}
	
	#pragma omp parallel for
	for( std::vector<std::vector<Matrixd>>::size_type idx = 0; idx < per_thread_C.size(); idx++ )
	{
		per_thread_C[idx].resize( NP );
		std::fill(per_thread_C[idx].begin(), per_thread_C[idx].end(), Matrixd::Zero());
	}
	
	#pragma omp parallel for
	for(int i=0; i<NP; ++i)
	{
		//get omp thread ID
		const int tid{ omp_get_thread_num() };
		assert( tid >= 0 );
		
		const Vectori base_grid =	((sim_state.mp.Position.col(i) - sim_state.SimulationBoxMin ) *width_inv - Vectord::Constant(0.5)).cast<int>();
		const Vectord d_from_basegrid = (sim_state.mp.Position.col(i) - sim_state.SimulationBoxMin )*width_inv - base_grid.cast<double>();
		
		//set basis function on omp thread ID
		basis_fs[tid].set(d_from_basegrid.cast<double>());
		
		const Vectori w = Vectori::Constant(3);
		const int Max = std::pow(3,SET::dim);
		for(int l=0; l<Max; ++l)
		{
			const Vectori node 			= MPM::flat2node(l, w);
			const Vectori base_node = base_grid + node;
			const Vectord distance_from_index = (node.cast<double>() - d_from_basegrid);

			per_thread_Velocity[tid].col(i) += basis_fs[tid].tensor(node)*sim_state.mg.Velocity.col( MPM::node2flat(base_node, sim_state.SimulationGridsnum)  ) ;

			per_thread_C[tid][i] += 4 * width_inv * basis_fs[tid].tensor(node) * sim_state.mg.Velocity.col(MPM::node2flat(base_node, sim_state.SimulationGridsnum) )*(distance_from_index).transpose();
		}
	}
	
	for( std::vector<MatrixXd>::size_type idx = 0; idx < per_thread_Velocity.size(); idx++ )
	{
		sim_state.mp.Velocity += per_thread_Velocity[idx];
	}
	
	#pragma omp parallel for
	for(int i=0; i<NP; ++i)
	{
		for( std::vector<std::vector<Matrixd>>::size_type idx = 0; idx < per_thread_Velocity.size(); idx++ )
		{
			sim_state.mp.C[i] += per_thread_C[idx][i];
		}
	}
	#endif
	
}

inline void MPM::Integrater::advection()
{
	
	const int NP=sim_state.mp.numpoints;
	
	#ifdef USEOPENMP
	#pragma omp parallel for
	#endif
	for(int i=0; i<NP; ++i)
	{
		sim_state.mp.Position.col(i) += sim_state.dt*sim_state.mp.Velocity.col(i);
	}
	

}


void MPM::Integrater::deformation_update()
{
	const int NP=sim_state.mp.numpoints;
	
	#ifdef USEOPENMP
	#pragma omp parallel for
	#endif
	for(int i=0; i<NP; ++i)
	{
	// MLS-MPM Deformation F-update
		const auto f = ( Matrixd::Identity() + sim_state.dt * sim_state.mp.C[i] );
		const Matrixd F = f * sim_state.mp.F[i];
		const Matrixd b = f * sim_state.mp.b[i] * f.transpose();
		
		// MLS-MPM F Plastic Correction
		auto [F_new, J_new, b_new] = sim_state.c->cal_F_placstic_correction(F, sim_state.mp.J[i], b);
		//std::cout << "J     " << sim_state.mp.J[i] << std::endl;
		//std::cout << "J_new " << J_new << std::endl;
		sim_state.mp.F[i] = std::move(F_new);
		sim_state.mp.J[i] = std::move(J_new);
		sim_state.mp.b[i] = std::move(b_new);
	}
	
	//sim_state.mp.show();
}

inline void MPM::Integrater::updateState()
{
	++sim_state.iterate_num;
}

inline void MPM::Integrater::stateCheck()
{
	const double timestep = sim_state.dt*sim_state.iterate_num;
	//std::cout << "time step:  " << timestep << std::endl;
	if(timestep > sim_state.simulation_end_time)
	{
		sim_state.simulation_stop = true;
	}
}



#endif /* mpm_integrator_hpp */
