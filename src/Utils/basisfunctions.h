//
//  basisfunctions.h
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2020/12/12.
//  Copyright © 2020年 Kentaro Nagasawa. All rights reserved.
//

#ifndef basisfunctions_h
#define basisfunctions_h

#include <memory>
#include <iostream>
#include <Eigen/Dense>

template <typename T, int dim> using Vector = Eigen::Matrix<T, dim, 1>;
template <typename T, int dim> using Matrix = Eigen::Matrix<T, dim, dim>;

template <typename T, int dim> using VectorX = Eigen::Matrix<T,Eigen::Dynamic,dim>;
template <typename T, int dim> using MatrixX = Eigen::Matrix<T, dim, Eigen::Dynamic>;
namespace MPM{

	template<int dim>
	inline auto flat2node(const int& i, const Vector<int, dim>& dim_width){
		if constexpr(dim == 2){
			//std::cout << "dim == 2" << std::endl;
			return Vector<int, 2>{i % dim_width.x(), (i / dim_width.x()) % dim_width.y()};
			
		}
		else if constexpr(dim == 3){
			//std::cout << "dim == 3" << std::endl;
			return Vector<int, 3>{i % dim_width.x(), (i / dim_width.x()) % dim_width.y(), i/(dim_width.x()*dim_width.y()) };
		}
		else{
			//static_assert(dim, "Cond must be false");
			std::cout << "dim == else" << std::endl;
		}
	}

	template<int dim>
	inline auto node2flat(const Vector<int, dim>& I, const Vector<int, dim>& dim_width){
		if constexpr(dim == 2){
			//std::cout << "dim == 2" << std::endl;
			double flat = I.y() * ( dim_width.x() ) + I.x();
			return flat;
		}
		else if constexpr(dim == 3){
			//std::cout << "dim == 3" << std::endl;
			double flat = I.z() * ( dim_width.y() ) * ( dim_width.x() ) + I.y() * ( dim_width.x() ) + I.x();
			return flat;
		}
		else{
			//static_assert(dim, "Cond must be false");
			std::cout << "dim == else" << std::endl;
		}
	}




}



using Array3d = Eigen::Array<double,3,1>;
using Array3u = Eigen::Array<unsigned,3,1>;
using Vector2d = Eigen::Matrix<double,2,1>;
using Vector3d = Eigen::Matrix<double,3,1>;
using Vector4d = Eigen::Matrix<double,4,1>;
using Vector3u = Eigen::Matrix<unsigned,3,1>;
using Vector6d = Eigen::Matrix<double,6,1>;
using VectorXd = Eigen::Matrix<double,Eigen::Dynamic,1>;
using VectorXb = Eigen::Matrix<bool,Eigen::Dynamic,1>;
using Matrix3d = Eigen::Matrix<double,3,3>;
using Matrix3Xd = Eigen::Matrix<double,3,Eigen::Dynamic>;
using Matrix4Xd = Eigen::Matrix<double,4,Eigen::Dynamic>;





#endif /* basisfunctions_h */

