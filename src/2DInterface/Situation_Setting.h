//
//  Situation_Setting.h
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2021/01/07.
//  Copyright © 2021 Kentaro Nagasawa. All rights reserved.
//

#ifndef Situation_Setting_h
#define Situation_Setting_h
#include "../SimulationCore/mpm_state.h"

namespace MPM{
//void Set_Material_Point_Randomly(const Vectord& center,const double& mat_width, const Vector3d& color,  const MPM::SimulationParameters& sim_params, MPM::SimulationState& sim_setting)


struct PointCube{
	Vectord min;
	Vectord max;
	double  rho;
	
	PointCube (){};
	PointCube (Vectord min, Vectord max, double rho)
		: min(min),
			max(max),
			rho(rho)
		{};
	void Set_MaterialPoints(MPM::SimulationState& sim_setting, const Vector3d& color);
};

struct SetMaterialPoint{
	std::vector<PointCube>		pointcubes;
};

void PointCube::Set_MaterialPoints(MPM::SimulationState& sim_setting, const Vector3d& color)
{
	Vectord small_val = Vectord::Constant(0.000001);
	Vectord padding = Vectord::Constant(sim_setting.SimulationGridsWidth/(sim_setting.materialpointsnum + 1));
	Vectori totalpoints = sim_setting.materialpointsnum*((max - min + small_val)/sim_setting.SimulationGridsWidth).cast<int>();
	double hl = 0.5*(sim_setting.SimulationGridsWidth/sim_setting.materialpointsnum);
	double volume = 8*hl*hl*hl;
	double mass = rho*volume;
	std::cout << mass << std::endl;
	int num = 1;
	for(int i=0;i<SET::dim;++i)
	{
		num *= totalpoints[i];
	}
	
	for(int i=0;i<num;++i)
	{
		auto node = MPM::flat2node<SET::dim>(i, totalpoints);
		const Vectord mp_pos = node.cast<double>()*(sim_setting.SimulationGridsWidth /sim_setting.materialpointsnum) + min + padding;
		sim_setting.mp.add_points(mp_pos, Vectord::Zero(SET::dim), mass, volume, color);
		//std::cout << i << " " << mp_pos << std::endl;
	}
	std::cout << "total points " << totalpoints << std::endl;
	std::cout << "num " << num << std::endl;
}


struct SituationSets{
	MPM::SetMaterialPoint			set_mat_point;
};


inline void Set_SimulationBoxSetting(const MPM::SimulationParameters& sim_params, MPM::SimulationState& sim_setting)
{
	//sim_setting.SimulationBoxMin = sim_params.box_min;
	//sim_setting.SimulationBoxMax = sim_params.box_max;
	
	sim_setting.SimulationGridsWidth		 = sim_params.d_width;
	sim_setting.SimulationGridsWidth_inv = 1.0/sim_params.d_width;
	
	const Vectord Boundary_padding = 5.0 * sim_setting.SimulationGridsWidth *Vectord::Constant(1.0);

	sim_setting.SimulationBoxMinBound = sim_params.box_min;
	sim_setting.SimulationBoxMaxBound = sim_params.box_max;

	sim_setting.SimulationBoxMin = sim_params.box_min - Boundary_padding;
	sim_setting.SimulationBoxMax = sim_params.box_max + Boundary_padding;

	
	sim_setting.SimulationBoxSize = sim_setting.SimulationBoxMax - sim_setting.SimulationBoxMin;
	
	sim_setting.SimulationGridsnum = (sim_setting.SimulationBoxSize / sim_setting.SimulationGridsWidth ).cast<int>();
	
	sim_setting.dt = sim_params.dt;
	sim_setting.simulation_end_time  = sim_params.end_time;
	
	ConstitutiveModels model_name;
	if ( sim_params.model_name == "Elastic"){
		model_name = ConstitutiveModels::Elastic;
	}
	else if ( sim_params.model_name == "SnowPlastic"){
		model_name = ConstitutiveModels::SnowPlastic;
	}
	else if ( sim_params.model_name == "HerschelBulkley"){
		model_name = ConstitutiveModels::HerschelBulkley;
	}
	else {
		std::cout << "invalied constitutive name " << std::endl; exit(0);
	}
	
	sim_setting.Set_ConstitutiveModel(model_name);

	std::cout << "**Set_SimulationBoxSetting**" << std::endl;
	std::cout << " ->SimulationBoxMin 		: " << sim_setting.SimulationBoxMin.x() << " " << sim_setting.SimulationBoxMin.y() << std::endl;
	std::cout << " ->SimulationBoxMax 		: " << sim_setting.SimulationBoxMax.x() << " " << sim_setting.SimulationBoxMax.y() << std::endl;
	std::cout << " ->SimulationBoxSize 		: " << sim_setting.SimulationBoxSize.x() << " " << sim_setting.SimulationBoxSize.y() << std::endl;
	std::cout << " ->SimulationGridsnum 	: " << sim_setting.SimulationGridsnum.x() << " " << sim_setting.SimulationGridsnum.y() << std::endl;
	std::cout << " ->SimulationGridsWidth : " << sim_setting.SimulationGridsWidth << std::endl;
	std::cout << " ->Simulation dt 				: " << sim_setting.dt << std::endl;
	std::cout << " ->Simulation end time	: " << sim_setting.simulation_end_time << std::endl;
	
}



inline void Set_Material_Point_Randomly(const Vectord& center,const double& mat_width, const Vector3d& color,  const MPM::SimulationParameters& sim_params, MPM::SimulationState& sim_setting)
{
	for (int i = 0; i < 300; i++) {
		//auto rand = Eigen::MatrixXd::Random(SET::dim, 1)*2.0 - Eigen::MatrixXd::One(SET::dim, 1);
		Vectord rand = Vectord::Random(SET::dim);
		const Vectord mp_pos_rand = rand*mat_width + center;
		sim_setting.mp.add_points(mp_pos_rand, Vectord::Zero(SET::dim),1.0, 1.0, color);
	}
}

/*
inline void Set_Material_Point_FromFile(MPM::SimulationState& sim_setting, const Vector3d& color)
{
	unsigned long num_inserted_points{ 0 };
	// Attempt to open the sdf file
	std::string pfile = "/Users/kn/Downloads/MPMTest/bunny_pt";
	Vector3d start_point = Vector3d::Zero();
	FILE *point_file_dat = 0;
	try
	{

	point_file_dat = fopen( pfile.c_str(), "rb" );
	//if( !point_file_dat ) return material_points;
	
 
	// Read int the number of point
	std::int32_t num_points;
	fread( &num_points, sizeof(std::int32_t), 1, point_file_dat );
	// Read int the volume of point
	double point_volume;
	fread( &point_volume, sizeof(double), 1, point_file_dat );
	// Read int hight of point
	double hl;
	fread( &hl, sizeof(double), 1, point_file_dat );
	//
	std::cout << "hl is " << hl << std::endl;

	for(int i = 0; i < num_points ; i++ ){
		double position_data[3];
		fread( position_data, sizeof(double), 3, point_file_dat );
		
		// Adopt -Z Foeward  :: -1*(position_data[2] + start_point.z()
		const Vector3d position(position_data[0] + start_point.x(), position_data[1] + start_point.y(), -1*(position_data[2] + start_point.z()) );
		//std::cout << position.x() << " " << position.y() << " " << position.z() << std::endl;
		sim_setting.mp.add_points(position, Vectord::Zero(SET::dim), color);
		
		num_inserted_points ++;
	}
	
	// Close the output file
	fclose( point_file_dat );
	//std::cout << "laod point file" << point_file.file_name_obj << "completed" << std::endl;
	}
	catch( const std::string& error )
	{
		std::cerr << error << std::endl;
		exit(0);
	}
		
		
}
*/

};

#endif /* Situation_Setting_h */
