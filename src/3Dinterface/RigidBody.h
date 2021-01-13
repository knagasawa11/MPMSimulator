//
//  RigidBody.h
//  MPMSimulator_3D
//
//  Created by Kentaro Nagasawa on 2021/01/07.
//

#ifndef RigidBody_h
#define RigidBody_h

#ifdef MPM3D
#include "../3DInterface/mpm_setting.h"
#endif

#include "../Utils/basisfunctions.h"

using Vectord = Vector<double, SET::dim>;
using Vectori = Vector<int, SET::dim>;
using Matrixd = Matrix<double, SET::dim>;
using Matrixi = Matrix<int, SET::dim>;

using MatrixXd = MatrixX<double, SET::dim>;
namespace MPM{
struct RigidBody{
	std::string m_filename;
	Vectori m_grid_dimensions;
	Vectord m_cell_delta;
	Vectord m_grid_start;
	Vectord m_start;
	//Matrix3Xsc m_sdf_normals;
	VectorXd m_sdf_vals;
	Matrix3Xd m_sdf_normals;
	
	RigidBody (){};
	RigidBody (std::string filename)
		: m_filename(filename)
		{};
	
	bool loadSDF(std::string sdf_file_name);
	bool loadSDF();
	bool getSDF(const Vectord& x, double& sdfValue, Vectord& n, double boundaryOffset);
	
	void computeNormals();
	bool computeSDF(const Vectord& x, double& sdfValue, Vectord& n, double boundaryOffset);

	// Interpolation functions
	double trilinearInterpolation( const double& v000, const double& v100, const double& v010, const double& v110, const double& v001, const double& v101, const double& v011, const double& v111, const Vector3d& bc ) const;
	double linearInterpolation( const double& v0, const double& v1, const double& alpha ) const;
	void trilinearInterpolation( const Vector3d& v000, const Vector3d& v100, const Vector3d& v010, const Vector3d& v110, const Vector3d& v001, const Vector3d& v101, const Vector3d& v011, const Vector3d& v111, const Vector3d& bc, Vector3d& result ) const;
	void linearInterpolation( const Vector3d& v0, const Vector3d& v1, const double& alpha, Vector3d& result ) const;
	//---------------------------

};

}

inline bool MPM::RigidBody::loadSDF()
{
	// Attempt to open the sdf file
	std::cout << "Attempt to open the sdf file : " << m_filename << std::endl;
	FILE *sdf_file = 0;
	sdf_file = fopen( m_filename.c_str(), "rb" );
	if( !sdf_file ) return false;
	
	// Read in the min point of the grid
	double grid_start[SET::dim];
	fread( grid_start, sizeof(double), SET::dim, sdf_file );
	// Read in the space between grid points
	double cell_delta[SET::dim];
	fread( cell_delta, sizeof(double), SET::dim, sdf_file );
	// Read in the number of grid points along each direction
	std::int32_t grid_dimension[SET::dim];
	fread( grid_dimension, sizeof(std::int32_t), SET::dim, sdf_file );
	
	std::cout << grid_start[0] << " " << grid_start[1] << " " << grid_start[2] << std::endl;
	std::cout << cell_delta[0] << " " << cell_delta[1] << " " << cell_delta[2] << std::endl;
	std::cout << grid_dimension[0] << " " << grid_dimension[1] << " " << grid_dimension[2] << std::endl;
	
	double flat_grid_dimension{1.0};
	for(int i=0; i<SET::dim; ++i)
	{
		flat_grid_dimension *= grid_dimension[i];
	}
	
	// Read in the grid values along each direction
	double* data = (double*)malloc( sizeof(double) * flat_grid_dimension);
	fread( data, sizeof(double), flat_grid_dimension, sdf_file );
	// Close the output file
	fclose( sdf_file );
	
	// setting
	for(int i=0; i<SET::dim; ++i){
		m_grid_dimensions(i) = grid_dimension[i];
		m_cell_delta(i) = cell_delta[i];
		m_grid_start(i) = grid_start[i];
	}
	m_grid_start(2) += 1.0;
																 
	m_sdf_vals.resize(flat_grid_dimension);
	
	for(int i=0; i<flat_grid_dimension; i++){
		//std::cout << i << " " << data[i] << std::endl;
		m_sdf_vals[i] = data[i];
	}
		
	m_sdf_normals.resize(Eigen::NoChange, flat_grid_dimension);
	computeNormals();

	//m_sdf_normals.resize(Eigen::NoChange, flat_grid_dimension);
	//computeNormals();

	//m_Center(0) = m_grid_start(0) + m_cell_delta(0) + m_grid_dimensions(0) * 0.5;
	//m_Center(1) = m_grid_start(1) + m_cell_delta(1) + m_grid_dimensions(1) * 0.5;
	//m_Center(2) = m_grid_start(2) + m_cell_delta(2) + m_grid_dimensions(2) * 0.5;
	//m_InitialCenter = m_Center;
	//m_Velocity = Vector3s::Zero();
	//m_Translation = Vector3s::Zero();
	//quat = Eigen::Quaternionf::Identity();
	
	free(data);
	
	std::cout << "  Loaded SDF: " << m_filename << std::endl;
	
	return true;
}


inline bool MPM::RigidBody::loadSDF(std::string sdf_file_name)
{
	// Attempt to open the sdf file
	std::cout << "Attempt to open the sdf file : " << sdf_file_name << std::endl;
	FILE *sdf_file = 0;
	sdf_file = fopen( sdf_file_name.c_str(), "rb" );
	if( !sdf_file ) return false;
	
	// Read in the min point of the grid
	double grid_start[SET::dim];
	fread( grid_start, sizeof(double), SET::dim, sdf_file );
	// Read in the space between grid points
	double cell_delta[SET::dim];
	fread( cell_delta, sizeof(double), SET::dim, sdf_file );
	// Read in the number of grid points along each direction
	std::int32_t grid_dimension[SET::dim];
	fread( grid_dimension, sizeof(std::int32_t), SET::dim, sdf_file );
	
	std::cout << grid_start[0] << " " << grid_start[1] << " " << grid_start[2] << std::endl;
	std::cout << cell_delta[0] << " " << cell_delta[1] << " " << cell_delta[2] << std::endl;
	std::cout << grid_dimension[0] << " " << grid_dimension[1] << " " << grid_dimension[2] << std::endl;
	
	double flat_grid_dimension{1.0};
	for(int i=0; i<SET::dim; ++i)
	{
		flat_grid_dimension *= grid_dimension[i];
	}
	
	// Read in the grid values along each direction
	double* data = (double*)malloc( sizeof(double) * flat_grid_dimension);
	fread( data, sizeof(double), flat_grid_dimension, sdf_file );
	// Close the output file
	fclose( sdf_file );
	
	// setting
	for(int i=0; i<SET::dim; ++i){
		m_grid_dimensions(i) = grid_dimension[i];
		m_cell_delta(i) = cell_delta[i];
		m_grid_start(i) = grid_start[i];
	}
	m_grid_start(2) += 1.0;
																 
	m_sdf_vals.resize(flat_grid_dimension);
	
	for(int i=0; i<flat_grid_dimension; i++){
		//std::cout << i << " " << data[i] << std::endl;
		m_sdf_vals[i] = data[i];
	}
		
	m_sdf_normals.resize(Eigen::NoChange, flat_grid_dimension);
	computeNormals();

	//m_sdf_normals.resize(Eigen::NoChange, flat_grid_dimension);
	//computeNormals();

	//m_Center(0) = m_grid_start(0) + m_cell_delta(0) + m_grid_dimensions(0) * 0.5;
	//m_Center(1) = m_grid_start(1) + m_cell_delta(1) + m_grid_dimensions(1) * 0.5;
	//m_Center(2) = m_grid_start(2) + m_cell_delta(2) + m_grid_dimensions(2) * 0.5;
	//m_InitialCenter = m_Center;
	//m_Velocity = Vector3s::Zero();
	//m_Translation = Vector3s::Zero();
	//quat = Eigen::Quaternionf::Identity();
	
	free(data);
	
	std::cout << "  Loaded SDF: " << sdf_file_name << std::endl;
	
	return true;
}

inline bool MPM::RigidBody::getSDF(const Vectord& x, double& sdfValue, Vectord& n, double boundaryOffset)
{
	Vector3d small_value;
	small_value << 0.0000001,0.0000001,0.0000001;
	const Vector3d& _x = x - (small_value);
	return computeSDF(_x, sdfValue, n, boundaryOffset);
}

void MPM::RigidBody::computeNormals()
{
	for( unsigned k = 0; k < m_grid_dimensions.z(); ++k )
	{
		for( unsigned j = 0; j < m_grid_dimensions.y(); ++j )
		{
			for( unsigned i = 0; i < m_grid_dimensions.x(); ++i )
			{
		unsigned im = (i>0) ? i-1 : i; unsigned ip = (i<m_grid_dimensions.x()-1) ? i+1 : i;
		unsigned jm = (j>0) ? j-1 : j; unsigned jp = (j<m_grid_dimensions.y()-1) ? j+1 : j;
		unsigned km = (k>0) ? k-1 : k; unsigned kp = (k<m_grid_dimensions.z()-1) ? k+1 : k;
		
		double dx = (ip-im) * m_cell_delta(0);
		double dy = (jp-jm) * m_cell_delta(1);
		double dz = (kp-km) * m_cell_delta(2);
		
		unsigned flat_idx_im = k * m_grid_dimensions.y() * m_grid_dimensions.x() + j * m_grid_dimensions.x() + im;
		double vim = m_sdf_vals[flat_idx_im];
		
		unsigned flat_idx_ip = k * m_grid_dimensions.y() * m_grid_dimensions.x() + j * m_grid_dimensions.x() + ip;
		double vip = m_sdf_vals[flat_idx_ip];
		
		unsigned flat_idx_jm = k * m_grid_dimensions.y() * m_grid_dimensions.x() + jm * m_grid_dimensions.x() + i;
		double vjm = m_sdf_vals[flat_idx_jm];
		
		unsigned flat_idx_jp = k * m_grid_dimensions.y() * m_grid_dimensions.x() + jp * m_grid_dimensions.x() + i;
		double vjp = m_sdf_vals[flat_idx_jp];
		
		unsigned flat_idx_km = km * m_grid_dimensions.y() * m_grid_dimensions.x() + j * m_grid_dimensions.x() + i;
		double vkm = m_sdf_vals[flat_idx_km];
		
		unsigned flat_idx_kp = kp * m_grid_dimensions.y() * m_grid_dimensions.x() + j * m_grid_dimensions.x() + i;
		double vkp = m_sdf_vals[flat_idx_kp];
		
		double nx = (vip - vim) / dx;
		double ny = (vjp - vjm) / dy;
		double nz = (vkp - vkm) / dz;
		
		double inv_length = 1.0 / sqrt(nx*nx+ny*ny+nz*nz);
		
		unsigned flat_idx = k * m_grid_dimensions.y() * m_grid_dimensions.x() + j * m_grid_dimensions.x() + i;
		m_sdf_normals(0, flat_idx) = nx * inv_length;
		m_sdf_normals(1, flat_idx) = ny * inv_length;
		m_sdf_normals(2, flat_idx) = nz * inv_length;
		}
		}
	}
}


inline bool MPM::RigidBody::computeSDF(const Vectord& x, double& v, Vectord& n, double boundaryOffset)
{
	Eigen::Matrix<int,3,1> indices = ( ( x - m_grid_start ).array() / m_cell_delta.array() ).unaryExpr( [](const double& y) { return floor(y); } ).cast<int>();
	
	if((indices(0) < 0) || (indices(1) < 0) || (indices(2) < 0) || (indices(0) >= (int)(m_grid_dimensions(0)-1)) || (indices(1) >= (int)(m_grid_dimensions(1)-1)) || (indices(2) >= (int)(m_grid_dimensions(2)-1)))
		return false;
	
	const Vector3d bc = ( x.array() - ( m_grid_start.array() + indices.cast<double>().array() * m_cell_delta.array() ) ).array() / m_cell_delta.array();
	
	assert( ( bc.array() >= 0.0 ).all() ); assert( ( bc.array() <= 1.0 ).all() );
	
	// Grab the value of the distance field at each grid point
	unsigned flat_idx_000 = indices.z() * m_grid_dimensions(0) * m_grid_dimensions(1) + indices.y() * m_grid_dimensions(0) + indices.x();
	unsigned flat_idx_100 = indices.z() * m_grid_dimensions(0) * m_grid_dimensions(1) + indices.y() * m_grid_dimensions(0) + indices.x()+1;
	unsigned flat_idx_010 = indices.z() * m_grid_dimensions(0) * m_grid_dimensions(1) + (indices.y()+1) * m_grid_dimensions(0) + indices.x();
	unsigned flat_idx_110 = indices.z() * m_grid_dimensions(0) * m_grid_dimensions(1) + (indices.y()+1) * m_grid_dimensions(0) + indices.x()+1;
	unsigned flat_idx_001 = (indices.z()+1) * m_grid_dimensions(0) * m_grid_dimensions(1) + indices.y() * m_grid_dimensions(0) + indices.x();
	unsigned flat_idx_101 = (indices.z()+1) * m_grid_dimensions(0) * m_grid_dimensions(1) + indices.y() * m_grid_dimensions(0) + indices.x()+1;
	unsigned flat_idx_011 = (indices.z()+1) * m_grid_dimensions(0) * m_grid_dimensions(1) + (indices.y()+1) * m_grid_dimensions(0) + indices.x();
	unsigned flat_idx_111 = (indices.z()+1) * m_grid_dimensions(0) * m_grid_dimensions(1) + (indices.y()+1) * m_grid_dimensions(0) + indices.x()+1;

	v = trilinearInterpolation( m_sdf_vals(flat_idx_000), m_sdf_vals(flat_idx_100), m_sdf_vals(flat_idx_010), m_sdf_vals(flat_idx_110),
			m_sdf_vals(flat_idx_001), m_sdf_vals(flat_idx_101), m_sdf_vals(flat_idx_011), m_sdf_vals(flat_idx_111), bc )
			- boundaryOffset;
	
	trilinearInterpolation( m_sdf_normals.col(flat_idx_000), m_sdf_normals.col(flat_idx_100), m_sdf_normals.col(flat_idx_010), m_sdf_normals.col(flat_idx_110),
			m_sdf_normals.col(flat_idx_001), m_sdf_normals.col(flat_idx_101), m_sdf_normals.col(flat_idx_011), m_sdf_normals.col(flat_idx_111), bc, n );
			
	n.normalize();
	
	// if the point of interpolated SDF is negative, return true.
	// it means there are something rigid objects at the point.
	return v < 0.0;
}

inline double MPM::RigidBody::trilinearInterpolation( const double& v000, const double& v100, const double& v010, const double& v110, const double& v001, const double& v101, const double& v011, const double& v111, const Vector3d& bc ) const
{
	assert( ( bc.array() >= 0.0 ).all() ); assert( ( bc.array() <= 1.0 ).all() );

	const double v00 = linearInterpolation( v000, v100, bc.x() );
	const double v10 = linearInterpolation( v010, v110, bc.x() );
	const double v01 = linearInterpolation( v001, v101, bc.x() );
	const double v11 = linearInterpolation( v011, v111, bc.x() );

	const double v0 = linearInterpolation( v00, v10, bc.y() );
	const double v1 = linearInterpolation( v01, v11, bc.y() );

	const double v = linearInterpolation( v0, v1, bc.z() );

	return v;
}

inline double MPM::RigidBody::linearInterpolation( const double& v0, const double& v1, const double& alpha ) const
{
	assert( alpha >= 0.0 ); assert( alpha <= 1.0 );
	return ( 1.0 - alpha ) * v0 + alpha * v1;
}

inline void MPM::RigidBody::trilinearInterpolation( const Vector3d& v000, const Vector3d& v100, const Vector3d& v010, const Vector3d& v110,const Vector3d& v001, const Vector3d& v101, const Vector3d& v011, const Vector3d& v111, const Vector3d& bc, Vector3d& result ) const
{
	assert( ( bc.array() >= 0.0 ).all() ); assert( ( bc.array() <= 1.0 ).all() );

	Vector3d v00; linearInterpolation( v000, v100, bc.x(), v00 );
	Vector3d v10; linearInterpolation( v010, v110, bc.x(), v10 );
	Vector3d v01; linearInterpolation( v001, v101, bc.x(), v01 );
	Vector3d v11; linearInterpolation( v011, v111, bc.x(), v11 );

	Vector3d v0; linearInterpolation( v00, v10, bc.y(), v0 );
	Vector3d v1; linearInterpolation( v01, v11, bc.y(), v1 );

	linearInterpolation( v0, v1, bc.z(), result );
}

inline void MPM::RigidBody::linearInterpolation( const Vector3d& v0, const Vector3d& v1, const double& alpha, Vector3d& result ) const
{
	assert( alpha >= 0.0 ); assert( alpha <= 1.0 );
	result = ( 1.0 - alpha ) * v0 + alpha * v1;
}


#endif /* RigidBody_h */
