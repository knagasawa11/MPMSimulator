//
//  mpm_grid.h
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2020/12/13.
//  Copyright © 2020年 Kentaro Nagasawa. All rights reserved.
//

#ifndef mpm_grid_h
#define mpm_grid_h
#include <vector>
#include "../Utils/basisfunctions.h"
#ifdef MPM2D
#include "../2DInterface/mpm_setting.h"
#endif
#ifdef MPM3D
#include "../3DInterface/mpm_setting.h"
#endif



using Vectord = Vector<double, SET::dim>;
using Vectori = Vector<int, SET::dim>;
using Matrixd = Matrix<double, SET::dim>;
using Matrixi = Matrix<int, SET::dim>;

using MatrixXd = MatrixX<double, SET::dim>;

/*
struct MaterialGrid {
	Eigen::Vector2d Position;
	Eigen::Vector2d Velocity;
	inline MaterialGrid( const Eigen::Vector2d& position, const Eigen::Vector2d& velocity);
	inline MaterialGrid (const MaterialGrid& rhs) : Position(rhs.Position), Velocity(rhs.Velocity)
	{
	}
};

inline MaterialGrid::MaterialGrid( const Eigen::Vector2d& position, const Eigen::Vector2d& velocity)
{
	Position = position;
	Velocity = velocity;
}
*/

namespace MPM {
	struct MaterialGrids {
		MatrixXd Position;
		MatrixXd Velocity;
		std::vector<double>	Mass;
		std::vector<int> 	ValiedGridNodeIndices;
		
		int numgrids{0};
		void add_grids(const Vectord& position, const Vectord& velocity);
	};
	
	inline void MaterialGrids::add_grids(const Vectord &position, const Vectord &velocity)
	{
		Position.conservativeResize( SET::dim, numgrids + 1 );
		Velocity.conservativeResize( SET::dim, numgrids + 1);
		
		Position.col(numgrids) = position;
		Velocity.col(numgrids) = velocity;
		Mass.push_back(0.0);
		++numgrids;
	}

}

#endif /* mpm_grid_h */
