//
//  main.cpp
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2020/12/12.
//  Copyright © 2020年 Kentaro Nagasawa. All rights reserved.
//

// 3 Dimension MPM Simulator
#define MPM3D

#include <iostream>
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#include <map>
#include <string>
#include <vector>
#include <time.h>
#include "omp.h"

#include "../SimulationCore/mpm_particle.h"
#include "../SimulationCore/mpm_state.h"
#include "../SimulationCore/mpm_integrator.hpp"
#include "../SimulationCore/constitutive_model.hpp"

#ifdef MPM2D
#include "../2DInterface/mpm_setting.h"
#endif

#ifdef MPM3D
#include "../3DInterface/mpm_setting.h"
#endif

#include "../3DInterface/RigidBody.h"
#include "../3DInterface/Situation_Setting.h"
#include "../3DInterface/opengl_view.hpp"


#include "../3DInterface/parser.h"
#include "../3DInterface/Json_Parse.hpp"
#include "../Utils/timer.h"
#include "../Utils/IO_util.h"


MPM::SimulationState 			sim_setting;
MPM::Integrater						integrater;
MPM::SituationSets				situation_set;
MPM::BoundaryCondition		boundary_conditions;
MPM::MaterialParameters		mat_params;
MPM::SimulationParameters sim_params;
IO::IOSettings						io_setting;
timer timer_main;

namespace MPM{
void ShowSimulationSituation();
}

int main(int argc, char *argv[]) {
	

	// Parse command line options
	if( !PARSE::parseCommandLineOptions( &argc, &argv, mat_params, io_setting ))
	{
		std::cerr << "Failed to Parse Command " << argv[1] << std::endl;
		exit(-1);
	}
	// Parse Json Setting file and Set parameters
	PARSE::readJson(io_setting.input_json_file, mat_params, sim_params, situation_set);
	
	io_setting.init(sim_params);
	MPM::Set_SimulationBoxSetting(sim_params, sim_setting);
	
	//Attach Material parametr to SimulationState's Constitutive model
	sim_setting.c->set_material_parameter(mat_params, sim_params);
	sim_setting.SimulationSettingInit();
	
	//Atouch Simulation state to Integrator
	integrater.set_integrater(sim_setting);
	
	//Set Material Points
	for (auto& e : situation_set.set_mat_point.pointfiles) {
		std::cout << "SetPointFromFile Filename : " << e.point_file_string << std::endl;
		e.Set_Material_Point_FromFile(sim_setting, Vector3d{1.0, 0.0, 0.0});
	}
	for (auto& e : situation_set.set_mat_point.pointcubes) {
		std::cout << "SetPointCube : " << std::endl;
		e.Set_MaterialPoints(sim_setting, Vector3d{1.0, 0.0, 0.0});
	}
	//situation_set.set_mat_point.pointcubes[0].Set_MaterialPoints(sim_setting, Vector3d{1.0, 0.0, 0.0});
	/**
	sim_setting.Set_MaterialPointsRandomly(Vector<double, SET::dim>{0.45,0.65,0.45}, 0.08, Vector3d{1.0, 0.0, 0.0}, sim_params);
	sim_setting.Set_MaterialPointsRandomly(Vector<double, SET::dim>{0.55,0.45,0.55}, 0.08, Vector3d{0.0, 1.0, 0.0}, sim_params);
	sim_setting.Set_MaterialPointsRandomly(Vector<double, SET::dim>{0.55,0.85,0.55}, 0.08, Vector3d{0.0, 0.0, 1.0}, sim_params);
	/**/
	
	
	MPM::ShowSimulationSituation();
	timer_main.Get_start();
	// Simulation loop start
	
	#ifdef USEGL
	std::cout << "--Use opengl view--" << std::endl;
	glutInit(&argc, argv);
	glutInitWindowSize(640, 480);
	glutCreateWindow(argv[0]);
	
	glutDisplayFunc(OPENGL::display);
	//glutMouseFunc(OPENGL::mouse);
	//glutMotionFunc(OPENGL::mousemove);
	glutReshapeFunc(OPENGL::resize);
	glutIdleFunc(OPENGL::idle);
	glutMainLoop();
	#else
	std::cout << "--Unuse opengl view--" << std::endl;
	while(1)
	{
		integrater.integrate();
		//glutPostRedisplay();
		if(sim_setting.iterate_num%100 == 0)
		{
			std::cout << sim_setting.iterate_num << std::endl;
			MPM::WRITE::show();
		}
		if(sim_setting.simulation_stop == true)
		{
			timer_main.Get_end();
			timer_main.Show_duration();
			exit(0);
		}

	}
	#endif
	
	
	timer_main.Get_end();
	timer_main.Show_duration();
	
	return 0;
}

void MPM::ShowSimulationSituation(){
	std::cout << "** Show Simulation Situation **" << std::endl;
	std::cout << " -> Total Material Points : " << sim_setting.mp.numpoints <<std::endl;
	std::cout << " -> Total Grid		 Points : " << sim_setting.mg.numgrids <<std::endl;
}



