//
//  parser.h
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2021/01/06.
//  Copyright © 2021 Kentaro Nagasawa. All rights reserved.
//

#ifndef parser_h
#define parser_h
#include <getopt.h>
#ifdef MPM2D
#include "../2DInterface/mpm_setting.h"
#endif

#ifdef MPM3D
#include "../3DInterface/mpm_setting.h"
#endif

#include "../Utils/IO_util.h"


//extern MPM::MaterialParameters mat_params;

namespace PARSE{

	template<class T>
	bool extractFromString( const std::string& str, T& res )
	{
		std::stringstream input_strm( str );
		input_strm >> res;
		return !input_strm.fail();
	}


	static bool parseCommandLineOptions( int* argc, char*** argv, MPM::MaterialParameters& mat_params, IO::IOSettings& io_settings  )
	{
		const struct option long_options[] =
		{
			{ "output_dir", required_argument, nullptr, 'o' },
			{ "frequency", required_argument, nullptr, 'f' },
			{ nullptr, 0, nullptr, 0 }
		};
		
		while( true )
		{
			int option_index = 0;
			const int c{ getopt_long( *argc, *argv, "e:h:s:o:j:", long_options, &option_index ) };
			if( c == -1 )
			{
				break;
			}
			switch( c )
			{
				case 'e':
				{
					std::cerr << "load mat_params.eta for Herschel Bulkley model" << optarg <<  std::endl;
					if( !extractFromString( optarg, mat_params.eta ) )
					{
						std::cerr << "Failed to read valuee for argument for -e." << std::endl;
						return false;
					}
					break;
				}
				case 'h':
				{
					std::cerr << "load mat_params.h for Herschel Bulkley model" << optarg <<  std::endl;
					if( !extractFromString( optarg, mat_params.h ) )
					{
						std::cerr << "Failed to read valuee for argument for -h. " << std::endl;
						return false;
					}
					break;
				}
				case 's':
				{
					std::cerr << "load mat_params.sigma_Y for Herschel Bulkley model" << optarg <<  std::endl;
					if( !extractFromString( optarg, mat_params.sigma_Y ) )
					{
						std::cerr << "Failed to read valuee for argument for -s." << std::endl;
						return false;
					}
					break;
				}
				case 'o':
				{
					std::cerr << "load output_dir_path" << optarg <<  std::endl;
					io_settings.output_dir_path = optarg;
					break;
				}
				case 'j':
				{
					std::cerr << "load input_json_file" << optarg <<  std::endl;
					io_settings.input_json_file = optarg;
					break;
				}
				default:
				{
					std::cerr << "This is a bug in the command line parser. Please file a report." << std::endl;
					return false;
				}
			}
		}
		
		return true;
	}

};

#endif /* parser_h */
