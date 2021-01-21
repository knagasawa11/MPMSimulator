//
//  mpm_setting.h
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2020/12/31.
//  Copyright © 2020 Kentaro Nagasawa. All rights reserved.
//

#ifndef mpm_setting_h
#define mpm_setting_h

namespace SET{
//Simulation Dimension
constexpr int 		dim 				= 2				;


//Initial Simulation Parameter
constexpr double 	box_width 	= 1.0			;
constexpr double 	box_height 	= 1.0			;
constexpr double	bound_pad		= 0.1			;
constexpr double 	dt 					= 0.00001	;
constexpr double 	d_width 		= 0.0125	;
constexpr double 	end_time 		= 1.0			;
constexpr int			num_points 	= 2				;
constexpr int 		frame				= 50			;


//Initial Material Parameter
namespace MAT_PARAM{
constexpr double E 						= 10000.0	;
constexpr double nu 					= 0.2			;
constexpr double hardening 		= 10.0		;
constexpr double mu 					= 10000.0	;
constexpr double kappa 				= 10000.0	;
constexpr double eta 					= 1.0			;
constexpr double h 						= 1.0			;
constexpr double sigma_Y 			= 0.01		;
};


};
namespace MPM{
struct MaterialParameters{
		double E;
		double nu;
		double hardening;
		double mu;
		double kappa;
		double eta;
		double h;
		double sigma_Y;
		MaterialParameters ()
			: E					(SET::MAT_PARAM::E)
			, nu 				(SET::MAT_PARAM::nu)
			, hardening	(SET::MAT_PARAM::hardening)
			, mu				(SET::MAT_PARAM::mu)
			, kappa			(SET::MAT_PARAM::kappa)
			, eta				(SET::MAT_PARAM::eta)
			, h					(SET::MAT_PARAM::h)
			, sigma_Y		(SET::MAT_PARAM::sigma_Y)
			{};
};

struct SimulationParameters{
		Vector2d box_min		= Vector2d::Constant(0.0);
		Vector2d box_max		= Vector2d::Constant(1.0);
		double 	box_width 	= 1.0			;
		double 	box_height 	= 1.0			;
		double	bound_pad		= 0.1			;
		double 	dt 					= 0.00001	;
		double	alpha 			= 0.95		;
		double 	d_width 		= 0.0125	;
		double 	end_time 		= 1.0			;
		int			num_points 	= 2				;
		int 		frame				= 50			;
	  std::string model_name;
		SimulationParameters ()
				: box_width			(SET::box_width)
				, box_height 		(SET::box_height)
				, bound_pad			(SET::bound_pad)
				, dt						(SET::dt)
				,	alpha					(0.95)
				, d_width				(SET::d_width)
				, end_time			(SET::end_time)
				, num_points		(SET::num_points)
				, frame					(SET::frame)
				{};
};


};


#endif /* mpm_setting_h */
