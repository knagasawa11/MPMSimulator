//
//  mpm_integrator.cpp
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2020/12/18.
//  Copyright © 2020年 Kentaro Nagasawa. All rights reserved.
//

#include <vector>
#include "mpm_integrator.hpp"

//#include <unsupported/Eigen/CXX11/Tensor>



/*
MPM::Basisfunctions::Basisfunctions()
{
	N.push_back( Vectord::Ones(SET::dim) );
	N.push_back( Vectord::Ones(SET::dim) );
	N.push_back( Vectord::Ones(SET::dim) );
	std::cout << "N's size is " << N.size() << std::endl;
}
*/
std::vector<Vectord> Basisfunction(const Vectord &d_bg)
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


void MPM::Integrater::set_integrater(MPM::SimulationState& _sim_setting)
{
	sim_state = _sim_setting;
	basis_f.init();
}

void MPM::Integrater::integrate()
{
	stateCheck(sim_state);
	
	particle2grid(sim_state);
	gridAdvance(sim_state);
	grid2particle(sim_state);
	advection(sim_state);
	deformation_update(sim_state);
	
	updateState(sim_state);
}

void MPM::Integrater::position_update()
{
sim_state.mp.Position.row(1).array() -= 0.01;
}


void MPM::Integrater::particle2grid(MPM::SimulationState& _sim_setting)
{
	const double width_inv = sim_state.SimulationGridsWidth_inv;
	sim_state.mg.Velocity = Eigen::MatrixXd::Zero(SET::dim,sim_state.mg.numgrids+1);
	std::fill(sim_state.mg.Mass.begin(), sim_state.mg.Mass.end(), 0.0);
	
	const int NP=sim_state.mp.numpoints;
	for(int i=0; i<NP; ++i)
	{
		//const Matrixd FP = sim_state.c->cal_particle_stress(sim_state.mp.F[i], sim_state.mp.F[i].determinant());
		const Matrixd FP = sim_state.c->cal_particle_stress(sim_state.mp.F[i], sim_state.mp.J[i],sim_state.mp.b[i]);

		
		constexpr double p_vol = 1.0;
		constexpr double particle_mass = 1.0;
		const Matrixd Stress = - (sim_state.dt * p_vol) * (4 * width_inv * width_inv * FP);
		const Matrixd Affine = Stress + particle_mass * sim_state.mp.C[i];
		
		//const auto N = MPM::BasisFunction(d_from_basegrid.cast<double>());
		
		const Vectori base_grid =	(sim_state.mp.Position.col(i)*width_inv - Vectord::Constant(0.5)).cast<int>();
		const Vectord d_from_basegrid = sim_state.mp.Position.col(i)*width_inv - base_grid.cast<double>();

		
		basis_f.set(d_from_basegrid.cast<double>());
		//auto N = basis_f.N;

		//basis_f.set(d_from_basegrid.cast<double>());
		const Vectord mv = sim_state.mp.Velocity.col(i) * particle_mass;
		
		///
		const Vectori w = Vectori::Constant(3);
		//const auto vec_basis = basis_f.vec();
		
		//VectorXd mg_mass = VectorXd::Zero(std::pow(3,SET::dim));

		//for(int l=0, Max = std::pow(3,SET::dim); l<Max; ++l){
		//	Vectori node 			= MPM::flat2node(l, w);
		//	Vectori base_node = base_grid + node;
		//	mg_mass(l) = MPM::node2flat(base_node, sim_state.SimulationGridsnum);
		//}

		
		//for(int i=0, Max = std::pow(3,SET::dim); i<Max; ++i)
		for(int l=0, Max = std::pow(3,SET::dim); l<Max; ++l)
		{
			const Vectori node 			= MPM::flat2node(l, w);
			const Vectori base_node = base_grid + node;
			const Vectord distance_from_index = (node.cast<double>() - d_from_basegrid) * sim_state.SimulationGridsWidth;

			sim_state.mg.Velocity.col( MPM::node2flat(base_node, sim_state.SimulationGridsnum) ) += basis_f.tensor(node) * ( mv + Affine * distance_from_index );
			sim_state.mg.Mass[ MPM::node2flat(base_node, sim_state.SimulationGridsnum) ] += particle_mass * basis_f.tensor(node);

			//sim_state.mg.Velocity.col( MPM::node2flat(base_node, sim_state.SimulationGridsnum) ) += vec_basis(l) * ( mv + Affine * distance_from_index );
			//sim_state.mg.Mass[ MPM::node2flat(base_node, sim_state.SimulationGridsnum) ] += particle_mass * vec_basis(l);
		}
		///
		
		 
		/*
		for(int l=0;l<3;++l){for(int m=0;m<3;m++){

			const Vectord distance_from_index = (Eigen::Vector2d{l, m} - d_from_basegrid) * sim_state.SimulationGridsWidth;
			const Vectori node{l,m};
			sim_state.mg.Velocity.col(sim_state.node2flat(base_grid.x() + l, base_grid.y() + m)) += N[l].x()*N[m].y() * ( mv + Affine * distance_from_index );
			sim_state.mg.Mass[sim_state.node2flat(base_grid.x() + l, base_grid.y() + m)] += particle_mass * N[l].x()*N[m].y();
			//std::cout << "N[l].x()" << N[l].x() << std::endl;
			//std::cout << "mv" << mv << std::endl;
			//std::cout << "Affine" << Affine << std::endl;
			//std::cout << "distance_from_index" << distance_from_index << std::endl;
		}}
		*/

	}
}

void MPM::Integrater::gridAdvance(MPM::SimulationState &sim_state)
{
	for(int i=0,NP=sim_state.mg.numgrids;i<NP;++i)
	{
		if(sim_state.mg.Mass[i] > 0)
		{
			sim_state.mg.Velocity.col(i) = sim_state.mg.Velocity.col(i)/sim_state.mg.Mass[i];
			Vectord gravity = Vectord::Zero();
			gravity.y() = -198;
			sim_state.mg.Velocity.col(i) += sim_state.dt * gravity;
			
			
			const double bd_bottom = sim_state.SimulationGridsWidth*5;
			if(sim_state.mg.Position.col(i).y() < bd_bottom )
			{
				//sim_state.mg.Velocity.col(i) = Vectord::Constant(0.0);
				sim_state.mg.Velocity.col(i).y() = std::max(0.0,sim_state.mg.Velocity.col(i).y());
			}

			if((sim_state.mg.Position.col(i).array() < bd_bottom).any() || (sim_state.mg.Position.col(i).array() > (1.0 - bd_bottom) ).any() )
			{
				if(sim_state.mg.Position.col(i).y() < bd_bottom )
				{
					continue;
				}
				sim_state.mg.Velocity.col(i) = Vectord::Constant(0.0);
			}
		}
		
	}
}

void MPM::Integrater::grid2particle(MPM::SimulationState &sim_state)
{
	const double width_inv = sim_state.SimulationGridsWidth_inv;

	sim_state.mp.Velocity = Eigen::MatrixXd::Zero(SET::dim,sim_state.mp.numpoints);
	std::fill(sim_state.mp.C.begin(), sim_state.mp.C.end(), Matrixd::Zero());
	
	for(int i=0,NP=sim_state.mp.numpoints; i<NP; ++i)
	{
		
		const Vectori base_grid =	(sim_state.mp.Position.col(i)*width_inv - Vectord::Constant(0.5)).cast<int>();
		const Vectord d_from_basegrid = sim_state.mp.Position.col(i)*width_inv - base_grid.cast<double>();
		
		//auto N = MPM::BasisFunction(d_from_basegrid.cast<double>());
		//const auto N = basis_f.N;
		
		basis_f.set(d_from_basegrid.cast<double>());
		
		const Vectori w = Vectori::Constant(3);
		//const auto vec_basis = basis_f.vec();
		
		//for(int i=0, Max = std::pow(3,SET::dim); i<Max; ++i)
		
		for(int l=0, Max = std::pow(3,SET::dim); l<Max; ++l)
		{
			const Vectori node 			= MPM::flat2node(l, w);
			const Vectori base_node = base_grid + node;
			const Vectord distance_from_index = (node.cast<double>() - d_from_basegrid);

			sim_state.mp.Velocity.col(i)+= basis_f.tensor(node)*sim_state.mg.Velocity.col( MPM::node2flat(base_node, sim_state.SimulationGridsnum)  ) ;
			
			sim_state.mp.C[i] += 4 * width_inv * basis_f.tensor(node) * sim_state.mg.Velocity.col(MPM::node2flat(base_node, sim_state.SimulationGridsnum) )*(distance_from_index).transpose();
		}
		
	/*
		for(int l=0;l<3;++l){
			for(int m=0;m<3;m++){
				Vectord distance_from_index = (Eigen::Vector2d{l, m} - d_from_basegrid);
				
				sim_state.mp.Velocity.col(i)+= N[l].x()*N[m].y()*sim_state.mg.Velocity.col(sim_state.node2flat(base_grid.x() + l, base_grid.y() + m)) ;
				
				sim_state.mp.C[i] += 4 * width_inv * (N[l].x() * N[m].y() * sim_state.mg.Velocity.col(sim_state.node2flat(base_grid.x() + l, base_grid.y() + m)))*(distance_from_index).transpose();
			}
		}
	 */
	 
		
	}
}

void MPM::Integrater::advection(MPM::SimulationState &sim_state)
{
	sim_state.mp.Position += sim_state.dt*sim_state.mp.Velocity;
}


void MPM::Integrater::deformation_update(MPM::SimulationState &sim_state)
{
	for(int i=0,NP=sim_state.mp.numpoints; i<NP; ++i)
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

void MPM::Integrater::updateState(MPM::SimulationState &sim_state)
{
	++sim_state.iterate_num;
}

void MPM::Integrater::stateCheck(MPM::SimulationState &sim_state)
{
	const double timestep = sim_state.dt*sim_state.iterate_num;
	//std::cout << "time step:  " << timestep << std::endl;
	if(timestep > sim_state.simulation_end_time)
	{
		sim_state.simulation_stop = true;
	}
}

