//
//  Json_Parse.hpp
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2021/01/06.
//  Copyright © 2021 Kentaro Nagasawa. All rights reserved.
//

#ifndef Json_Parse_hpp
#define Json_Parse_hpp
#include <iostream>
#include <fstream>
#include "nlohmann/json.hpp"
#ifdef MPM2D
#include "../2DInterface/mpm_setting.h"
#endif

#ifdef MPM3D
#include "../3DInterface/mpm_setting.h"
#endif




namespace PARSE{

template <typename T>
T CheckJson(nlohmann::json j, std::string str1){
	if(j[str1].is_null()){
		std::cout << "Invalied parse: " << str1 << std::endl;
	}
	return j[str1].get<T>();
}
template <typename T>
T CheckJson(nlohmann::json j, std::string str1, std::string str2){
	if(j[str1][str2].is_null()){
		std::cout << "Invalied parse: " << str1 << std::endl;
	}
	return j[str1][str2].get<T>();
}

template <>
Vector2d CheckJson<Vector2d>(nlohmann::json j, std::string str1, std::string str2){
	if(j[str1][str2].is_null()){
		std::cout << "Invalied parse: " << str1 << std::endl;
	}
	return Vector2d{ j[str1][str2][0].get<double>(), j[str1][str2][1].get<double>()};
}

template <>
Vector2d CheckJson<Vector2d>(nlohmann::json j, std::string str1){
	if(j[str1].is_null()){
		std::cout << "Invalied parse: " << str1 << std::endl;
	}
	return Vector2d{ j[str1][0].get<double>(), j[str1][1].get<double>() };
}




inline void readJson(const std::string json_file, MPM::MaterialParameters& mat_params, MPM::SimulationParameters& sim_params, MPM::SituationSets& situation_set)
{
	std::cout << ":: Read Json ::" << std::endl;
	std::ifstream reading(json_file);
	nlohmann::json j;
	reading >> j;

	
	// Read Simulation Setting from Json
	std::cout << "SimSetting" << std::endl;
	
	std::cout << "-> box_min	: " << 	j["SimSetting"]["box_min"] 		<< std::endl;
	sim_params.box_min = CheckJson<Vector2d>(j,"SimSetting","box_min");
	
	std::cout << "-> box_max	: " << 	j["SimSetting"]["box_max"] 		<< std::endl;
	sim_params.box_max = CheckJson<Vector2d>(j,"SimSetting","box_max");
	
	std::cout << "-> box_width	: " << 	j["SimSetting"]["box_width"] 		<< std::endl;
	sim_params.box_width = CheckJson<double>(j,"SimSetting","box_width");
	
	std::cout << "-> bound_pad	: " << 	j["SimSetting"]["bound_pad"] 		<< std::endl;
	sim_params.bound_pad = CheckJson<double>(j,"SimSetting","bound_pad");
	
	std::cout << "-> dt					: " << 	j["SimSetting"]["dt"] 					<< std::endl;
	sim_params.dt = CheckJson<double>(j,"SimSetting","dt");
	
	std::cout << "-> alpha			: " << 	j["SimSetting"]["alpha"] 					<< std::endl;
	sim_params.alpha = CheckJson<double>(j,"SimSetting","alpha");
	
	std::cout << "-> d_width		: " << 	j["SimSetting"]["d_width"] 			<< std::endl;
	sim_params.d_width = CheckJson<double>(j,"SimSetting","d_width");

	std::cout << "-> end_time		: " << 	j["SimSetting"]["end_time"] 		<< std::endl;
	sim_params.end_time = CheckJson<double>(j,"SimSetting","end_time");

	std::cout << "-> num_points	: " << 	j["SimSetting"]["num_points"] 	<< std::endl;
	sim_params.num_points = CheckJson<int>(j,"SimSetting","num_points");

	std::cout << "-> frame			: " << 	j["SimSetting"]["frame"]				<< std::endl;
	sim_params.frame = CheckJson<int>(j,"SimSetting","frame");
	
	
	// Read Material Setting from Json
	std::cout << "MaterialCondition" << std::endl;
	std::cout << "-> ConstitutiveModel	: " << 	j["MaterialCondition"]["ConstitutiveModel"] 		<< std::endl;
	sim_params.model_name = CheckJson<std::string>(j,"MaterialCondition","ConstitutiveModel");
	
	std::cout << "-> E					: " << 	j["MaterialCondition"]["E"] 		<< std::endl;
	mat_params.E = CheckJson<double>(j,"MaterialCondition","E");
	
	std::cout << "-> nu					: " << 	j["MaterialCondition"]["nu"] 		<< std::endl;
	mat_params.nu = CheckJson<double>(j,"MaterialCondition","nu");
	
	std::cout << "-> hardening	: " << 	j["MaterialCondition"]["hardening"] 					<< std::endl;
	mat_params.hardening = CheckJson<double>(j,"MaterialCondition","hardening");
	
	std::cout << "-> mu					: " << 	j["MaterialCondition"]["mu"] 			<< std::endl;
	mat_params.mu = CheckJson<double>(j,"MaterialCondition","mu");
	
	std::cout << "-> kappa			: " << 	j["MaterialCondition"]["kappa"] 		<< std::endl;
	mat_params.kappa = CheckJson<double>(j,"MaterialCondition","kappa");
	
	std::cout << "-> eta				: " << 	j["MaterialCondition"]["eta"] 	<< std::endl;
	mat_params.eta = CheckJson<double>(j,"MaterialCondition","eta");
	
	std::cout << "-> h					: " << 	j["MaterialCondition"]["h"]				<< std::endl;
	mat_params.h = CheckJson<double>(j,"MaterialCondition","h");
	
	std::cout << "-> sigma_Y		: " << 	j["MaterialCondition"]["sigma_Y"]				<< std::endl;
	mat_params.sigma_Y = CheckJson<double>(j,"MaterialCondition","sigma_Y");
	
	
	std::cout << "SetPointCube" << std::endl;
	for (auto& element : j["SetPointCube"]) {
		Vector2d	min = CheckJson<Vector2d>(element,"min");
		Vector2d	max = CheckJson<Vector2d>(element,"max");
		double		rho = CheckJson<double>(element,"rho");
		situation_set.set_mat_point.pointcubes.push_back(MPM::PointCube{min, max, rho});
		//situation_set.set_mat_point.pointfiles.push_back(MPM::PointFile{filename});
	}
	
	/*
	std::cout << "Boundary	: "<< j["Boundary"].size()<< std::endl;
	for (auto& element : j["Boundary"]) {
		
		std::cout << element["position"] << '\n';
		std::cout << element["condition"] << '\n';
		
	}
	*/
	
}

}
#endif /* Json_Parse_hpp */
