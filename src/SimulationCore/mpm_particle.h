//
//  mpm_particle.h
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2020/12/12.
//  Copyright © 2020年 Kentaro Nagasawa. All rights reserved.
//

#ifndef mpm_particle_h
#define mpm_particle_h

#include "../Utils/basisfunctions.h"

#ifdef MPM2D
#include "../2DInterface/mpm_setting.h"
#endif

#ifdef MPM3D
#include "../3DInterface/mpm_setting.h"
#endif



//constexpr int  dim = 2;
using Vectord = Vector<double, SET::dim>;
using Vectori = Vector<int, SET::dim>;
using Matrixd = Matrix<double, SET::dim>;
using Matrixi = Matrix<int, SET::dim>;

using MatrixXd = MatrixX<double, SET::dim>;

/*
struct MaterialPoint {
	//Vector<double, dim> Position;
	Vectord Position;
	Vectord Velocity;
	inline MaterialPoint( const Vectord& position, const Vectord& velocity);
	inline MaterialPoint (const MaterialPoint& rhs) : Position(rhs.Position), Velocity(rhs.Velocity)
	{
	}
};

inline MaterialPoint::MaterialPoint( const Vectord& position, const Vectord& velocity)
{
	Position = position;
	Velocity = velocity;
}
*/

namespace MPM {
	struct MaterialPoints {
		//Vector<double, dim> Position;
		MatrixXd Position;
		MatrixXd Velocity;
		std::vector<Matrixd> F;
		std::vector<Matrixd> b;
		std::vector<double> J;
		std::vector<double> m;
		std::vector<double> vol;
		std::vector<Matrixd> C;
		std::vector<Vector3d> color;
		
		int numpoints{0};
		
		void add_points(const Vectord &position, const Vectord &velocity, const double& mass, const double& volume, const Vector3d& c);
		void show();
	};
	
	inline void MaterialPoints::add_points(const Vectord &position, const Vectord &velocity, const double& mass, const double& volume, const Vector3d& c)
	{
		Position.conservativeResize( SET::dim, numpoints + 1 );
		Velocity.conservativeResize( SET::dim, numpoints + 1);
		Position.col(numpoints) = position;
		Velocity.col(numpoints) = velocity;
		J.push_back(1.0);
		m.push_back(mass);
		vol.push_back(volume);
		F.push_back(Matrixd::Identity());
		b.push_back(Matrixd::Identity());
		C.push_back(Matrixd::Zero());
		color.push_back(c);
		++numpoints;
	}
	
	inline void MaterialPoints::show()
	{
		for(int i=0;i<numpoints;++i)
		{
			std::cout << "Point Nuber is " << i << std::endl;
			std::cout << "Position" << std::endl;
			std::cout << Position.col(i).x() << " " << Position.col(i).y() << std::endl;
			std::cout << "Velocity" << std::endl;
			std::cout << Velocity.col(i).x() << " " << Velocity.col(i).y() << std::endl;
			std::cout << "F" << std::endl;
			std::cout << F[i](0) << " " << F[i](1) << " " << F[i](2) << " " << F[i](3) << std::endl;
			std::cout << "C" << std::endl;
			std::cout << C[i](0) << " " << C[i](1) << " " << C[i](2) << " " << C[i](3) << std::endl;
			std::cout << "J" << std::endl;
			std::cout << J[i] << std::endl;
		}
	}

	
}

#endif /* mpm_particle_h */
