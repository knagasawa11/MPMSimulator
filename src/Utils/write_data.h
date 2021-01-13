//
//  write_data.h
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2021/01/05.
//  Copyright © 2021 Kentaro Nagasawa. All rights reserved.
//

#ifndef write_data_h
#define write_data_h
#include <filesystem>

#include "../SimulationCore/mpm_state.h"
#include "../Utils/IO_util.h"


extern MPM::SimulationState 	sim_setting;
extern IO::IOSettings					io_setting;

namespace MPM{
namespace WRITE{

static std::string generateFileName( const std::string prefix, const std::string extension, const std::string path, const unsigned int frame)
{
	std::stringstream ss;
	if( !path.empty() )
	{
		ss << path << "/";
	}
//		ss << prefix << "_" << std::setfill('0') << std::setw( int(g_save_number_width) ) << g_output_frame << "." << extension;
			ss << prefix << "_" << std::setfill('0') << std::setw( int(io_setting.output_num_width) ) << frame << "." << extension;
	return ss.str();
}



inline void show(){
	//std::cout << generateFileName( "config", "dat", io_setting.output_dir_path, io_setting.output_frame) << std::endl;
}

	inline void save_dat_file()
{
	// Write dat file for particle skinner
	const std::string output_file_name_dat{ generateFileName( "config", "dat", io_setting.output_dir_path, io_setting.output_frame) };
	std::cout << "Write dat file at " << output_file_name_dat << std::endl;
	/**/
	FILE *dat = 0;
	try
	{
		dat = fopen(output_file_name_dat.c_str(), "wb");
		
		int numpoints = sim_setting.mp.numpoints;
		int id = 1;
		float radi = 1.0;
		
		float data[3];
		float _data;
		std::int32_t _idata;
		
		//output numpoint
		_idata = numpoints;
		fwrite( &_idata, sizeof(std::int32_t), 1, dat );
		 
		if constexpr(SET::dim == 2){
			
			std::cout << "dim 2" << std::endl;
			//output point q
			for(int i = 0;i <numpoints;i++ ){
				data[0] = sim_setting.mp.Position.col(i).x();data[1] = sim_setting.mp.Position.col(i).y();data[2] = 0.0;
				fwrite( data, sizeof(float), 3, dat );
			}
			//output point radi
			for(int i = 0;i <numpoints;i++ ){
				_data = radi;
				fwrite( &_data, sizeof(float), 1, dat );
			}
			//output point velocity
			for(int i = 0;i <numpoints;i++ ){
				data[0] = sim_setting.mp.Velocity.col(i).x();data[1] = sim_setting.mp.Velocity.col(i).y();data[2] = 0.0;
				fwrite( data, sizeof(float), 3, dat );
			}
			//output point id
			for(int i = 0;i <numpoints;i++ ){
				_idata = id;
				fwrite( &_idata, sizeof(std::int32_t), 1, dat );
			}

		} else if constexpr(SET::dim == 3) {

			std::cout << "dim 3" << std::endl;
			//output point q
			for(int i = 0;i <numpoints;i++ ){
				data[0] = sim_setting.mp.Position.col(i).x();data[1] = sim_setting.mp.Position.col(i).y();data[2] = sim_setting.mp.Position.col(i).z();
				fwrite( data, sizeof(float), 3, dat );
			}
			//output point radi
			for(int i = 0;i <numpoints;i++ ){
				_data = radi;
				fwrite( &_data, sizeof(float), 1, dat );
			}
			//output point velocity
			for(int i = 0;i <numpoints;i++ ){
				data[0] = sim_setting.mp.Velocity.col(i).x();data[1] = sim_setting.mp.Velocity.col(i).y();data[2] = sim_setting.mp.Velocity.col(i).z();
				fwrite( data, sizeof(float), 3, dat );
			}
			//output point id
			for(int i = 0;i <numpoints;i++ ){
				_idata = id;
				fwrite( &_idata, sizeof(std::int32_t), 1, dat );
			}

		}
		fclose( dat );

	}
	catch( const std::string& error )
	{
		std::cerr << error << std::endl;
		exit(0);
	}
	/**/
}

}
}

#endif /* write_data_h */
