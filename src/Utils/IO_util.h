//
//  IO_util.h
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2021/01/06.
//  Copyright © 2021 Kentaro Nagasawa. All rights reserved.
//

#ifndef IO_util_h
#define IO_util_h
#include <iostream>

#ifdef MPM2D
#include "../2DInterface/mpm_setting.h"
#endif

#ifdef MPM3D
#include "../3DInterface/mpm_setting.h"
#endif

namespace IO
{

inline unsigned computeNumDigits( unsigned n )
{
	if( n == 0 ) { return 1; }
	unsigned num_digits{ 0 };
	while( n != 0 )
	{
		n /= 10;
		++num_digits;
	}
	return num_digits;
}

struct IOSettings
{
	std::string 	output_dir_path	;
	std::string 	input_json_file	;
	
	unsigned int 	output_frame		;
	int						output_num_width;
	int 				output_per_frame;
	void init(const MPM::SimulationParameters& sim_params);
	
	IOSettings ()
		: output_dir_path		(output_dir_path)
		, input_json_file		(input_json_file)
		, output_frame 			(0)
		, output_num_width	(1)
		, output_per_frame	(1)
		{};
};

inline void IOSettings::init(const MPM::SimulationParameters& sim_params){

	output_per_frame = (int)std::round(1.0/sim_params.dt)/sim_params.frame;
	output_num_width = computeNumDigits( 1 + unsigned( ceil(  sim_params.end_time / double( sim_params.dt ) ) ) / output_per_frame );

	std::cout << "input_json_file		:" << input_json_file << std::endl;
	std::cout << "output_dir_path		:" << output_dir_path << std::endl;
	std::cout << "output_frame			:" << output_frame << std::endl;
	std::cout << "output_num_width	:" << output_num_width << std::endl;
	std::cout << "output_per_frame	:" << output_per_frame << std::endl;
}

};

#endif /* IO_util_h */
