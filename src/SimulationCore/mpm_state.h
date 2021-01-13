//
//  mpm_state.h
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2020/12/13.
//  Copyright © 2020年 Kentaro Nagasawa. All rights reserved.
//

#ifndef mpm_state_h
#define mpm_state_h
#include <vector>

#include "mpm_grid.h"
#include "mpm_particle.h"
#include "constitutive_model.hpp"

#ifdef MPM2D
#include "../2DInterface/mpm_setting.h"
#endif

#ifdef MPM3D
#include "../3DInterface/mpm_setting.h"
#endif

#include "../Utils/basisfunctions.h"




namespace MPM{
	using Vectord = Vector<double, SET::dim>;
	using Vectori = Vector<int, SET::dim>;
	using Matrixd = Matrix<double, SET::dim>;
	using Matrixi = Matrix<int, SET::dim>;

	struct SimulationState {
		
		SimulationState() = default;
		
		static inline Vectord SimulationBoxMin;
		static inline Vectord SimulationBoxMax;
		static inline Vectord SimulationBoxMinBound;
		static inline Vectord SimulationBoxMaxBound;

		
		static inline Vectord SimulationBoxSize;
		static inline Vectori SimulationGridsnum;
		static inline double SimulationGridsWidth;
		static inline double SimulationGridsWidth_inv;
		
		//static inline std::vector<MaterialGrid> SimulationGrids;
		//static inline std::vector<MaterialPoint> MaterialPoints;
		static inline double dt;
		static inline int iterate_num{0};
		

		static inline MPM::MaterialPoints mp;
		static inline MPM::MaterialGrids  mg;
		static inline Elastic el;
		static inline ConstitutiveModel *c;

		int materialpointsnum{2};
		static inline double simulation_end_time;
		static inline bool simulation_stop;
		
		void Set_ConstitutiveModel(const ConstitutiveModels& model_name);
		void Set_SimulationSetting(const double & _dt, const double& _end_time);
		//void Set_SimulationBoxSize(const double& width, const double& height);
		void Set_SimulationBoxSize(const double& width);
		//void Set_SimulationGridsnum(const int& width,const int& height);
		void Set_SimulationGridsWidth (const double& width );
		void Set_SimulationGrids();
		//void Set_MaterialPoints(const double& minx,const double& Maxx,const double& miny,const double& Maxy);
		void Set_MaterialPoints(const Vectord& center,const double& mat_width, const MPM::SimulationParameters& sim_params);
		void Set_MaterialPointsRandomly(const Vectord& center,const double& mat_width, const Vector3d& color,  const MPM::SimulationParameters& sim_params);
		void SimulationSettingInit();
		void Show_var() const;
		
		static inline auto flat2node = [&a=SimulationGridsnum](int b){
			return Eigen::Vector2i{b % a.x(), (b / a.x()) % a.y()};
		};
		static inline auto node2flat = [&a=SimulationGridsnum](int i, int j){
			return int{i + (j * a.x())};
		};
		
	};



	inline void SimulationState::Set_ConstitutiveModel(const ConstitutiveModels &model_name)
	{
		switch (model_name) {
			case ConstitutiveModels::Elastic: c = new Elastic; break;
			case ConstitutiveModels::SnowPlastic: c = new SnowPlastic; break;
			case ConstitutiveModels::HerschelBulkley: c = new HerschelBulkley; break;
			default: std::cout << "invalied constitutive name " << std::endl; exit(0); break;
		}
	}

	inline void SimulationState::Set_SimulationSetting(const double& _dt, const double& _end_time)
	{
		dt = _dt;
		simulation_end_time = _end_time;
	}

	/*
	inline void SimulationState::Set_SimulationBoxSize( const double& width, const double& height)
	{
		SimulationBoxSize = Eigen::Vector2d{width, height};
	}
	 */

	inline void SimulationState::Set_SimulationBoxSize( const double& width)
	{
		SimulationBoxSize = width*Vectord::Ones(SET::dim);
	}

	/*
	inline void SimulationState::Set_SimulationGridsnum( const int& width, const int& height)
	{
		SimulationGridsnum = Eigen::Vector2i{width, height};
	}
	 */

	inline void SimulationState::Set_SimulationGridsWidth ( const double& width )
	{
		assert(width > 0);
		//assert(width < SimulationBoxSize.x() || width < SimulationBoxSize.y());
		
		SimulationGridsWidth = width;
		SimulationGridsWidth_inv = 1.0/width;
		//SimulationGridsnum = Vector<int,SET::dim>{SimulationBoxSize.x() / width, SimulationBoxSize.y() / width};
		SimulationGridsnum = (SimulationBoxSize / width ).cast<int>();
	}

	/*
	inline void SimulationState::Set_MaterialPoints(const double& minx,const double& Maxx,const double& miny,const double& Maxy)
	{
		const int totalpointsx = materialpointsnum*(int)((Maxx - minx)/SimulationGridsWidth);
		const int totalpointsy = materialpointsnum*(int)((Maxy - miny)/SimulationGridsWidth);
		for(int i = 0;i < totalpointsx; ++i)
		{
			for(int j = 0;j < totalpointsy; ++j){
				MaterialPoints.push_back(MaterialPoint(Eigen::Vector2d{minx + i*(SimulationGridsWidth/materialpointsnum), miny + j*(SimulationGridsWidth/materialpointsnum)}, Eigen::Vector2d{0.0, 0.0}));
				
				mp.add_points(Eigen::Vector2d{minx + i*(SimulationGridsWidth/materialpointsnum),miny + j*(SimulationGridsWidth/materialpointsnum)}, Eigen::Vector2d{0.0, -1.1});
			}
		}
	}
 */

	inline void SimulationState::Set_MaterialPoints(const Vectord& center,const double& mat_width, const MPM::SimulationParameters& sim_params)
	{
		Vectord min = center - mat_width*Vectord::Ones(SET::dim);
		Vectord max = center + mat_width*Vectord::Ones(SET::dim);
		Vectori totalpoints = materialpointsnum*((max - min)/SimulationGridsWidth).cast<int>();
		
		int num = 1;
		for(int i=0;i<SET::dim;++i)
		{
			num *= totalpoints[i];
		}
		
		for(int i=0;i<num;++i)
		{
			auto node = MPM::flat2node<SET::dim>(i, totalpoints);
			const Vectord mp_pos = node.cast<double>()*(sim_params.d_width /sim_params.num_points) + min;
			//mp.add_points(mp_pos,Vectord::Zero(SET::dim));
			
			//mp.Velocity.col(mp.numpoints - 1).x() = -1.0;
			//mp.add_points(Eigen::Vector2d{minx + i*(SimulationGridsWidth/materialpointsnum),miny + j*(SimulationGridsWidth/materialpointsnum)}, Eigen::Vector2d{0.0, -1.1});

		}
	}

	inline void SimulationState::Set_MaterialPointsRandomly(const Vectord& center,const double& mat_width, const Vector3d& color,  const MPM::SimulationParameters& sim_params)
	{
		for (int i = 0; i < 300; i++) {
			//auto rand = Eigen::MatrixXd::Random(SET::dim, 1)*2.0 - Eigen::MatrixXd::One(SET::dim, 1);
			Vectord rand = Vectord::Random(SET::dim);
			const Vectord mp_pos_rand = rand*mat_width + center;
			mp.add_points(mp_pos_rand, Vectord::Zero(SET::dim), color);
		}
	}


	inline void SimulationState::Set_SimulationGrids()
	{
		//const int total_grids_num = SimulationGridsnum.x()*SimulationGridsnum.y();
		int total_grids_num = 1;
		for(int i=0;i<SET::dim;++i)
		{
			total_grids_num *= SimulationGridsnum[i];
		}
		for (int i=0; i<total_grids_num; ++i)
		{
			//auto I = flat2node(i);
			auto I = MPM::flat2node<SET::dim>(i, SimulationGridsnum);
			const Vectord mg_pos = I.cast<double>()*SimulationGridsWidth + SimulationBoxMin;
			//auto p = Eigen::Vector2d{I.x(), i};
			//SimulationGrids.push_back(MaterialGrid(Eigen::Vector2d{I.x()*SimulationGridsWidth, I.y()*SimulationGridsWidth}, Eigen::Vector2d{0.0, 0.0}));
			//SimulationGrids.push_back(MaterialGrid(mg_pos, Vectord::Zero(SET::dim)));
			//mg.add_grids(Eigen::Vector2d{I.x()*SimulationGridsWidth, I.y()*SimulationGridsWidth}, Eigen::Vector2d{0.0, 0.0});
			mg.add_grids(mg_pos, Vectord::Zero(SET::dim));
		}
	}

	inline void SimulationState::SimulationSettingInit()
	{
		//Show_var();
		Set_SimulationGrids();
		simulation_stop = false;
		//Show_var();
	}
	/*
	inline void SimulationState::Show_var() const
	{
		std::cout << "SimulationBoxSize: " <<  SimulationBoxSize.x() << " " << SimulationBoxSize.y() << std::endl;
		std::cout << "SimulationGridsnum: " <<  SimulationGridsnum.x() << " " << SimulationGridsnum.y() << std::endl;
		std::cout << "SimulationGrids" << std::endl;
		for(auto g: SimulationGrids)
		{
			std::cout << g.Position.x() << " " << g.Position.y() << std::endl;
			std::cout << g.Velocity.x() << " " << g.Velocity.y() << std::endl;
		}
		std::cout << "MaterialPoints" << std::endl;
		for(auto g: MaterialPoints)
		{
			std::cout << g.Position.x() << " " << g.Position.y() << std::endl;
			std::cout << g.Velocity.x() << " " << g.Velocity.y() << std::endl;
		}
	}
	 */

}

#endif /* mpm_state_h */
