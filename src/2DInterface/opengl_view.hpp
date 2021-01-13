//
//  opengl_view.hpp
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2020/12/31.
//  Copyright © 2020 Kentaro Nagasawa. All rights reserved.
//

#ifndef opengl_view_hpp
#define opengl_view_hpp

#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#include <iostream>
#include <time.h>

#ifdef MPM2D
#include "../2DInterface/mpm_setting.h"
#endif


#include "../SimulationCore/mpm_state.h"
#include "../SimulationCore/mpm_integrator.hpp"

#include "../Utils/basisfunctions.h"

#include "../Utils/timer.h"
#include "../Utils/write_data.h"

extern MPM::SimulationState 	sim_setting;
extern MPM::Integrater				integrater;
extern IO::IOSettings					io_setting;
extern timer timer_main;



namespace OPENGL{

	inline double trans(double val)
	{
		return ((3*val/(sim_setting.SimulationBoxSize.x())) - 1.0);
	}

	template<int dim>
	inline void draw_dot(const Vectord& vec,const double& ps )
	{
		
		glPointSize(ps);
		glLineWidth(2.5);
		//glColor3f(1.0, 0.0, 0.0);
		glBegin(GL_POINTS);
		if constexpr(dim == 2){
		glVertex3f(trans(vec.x()),trans(vec.y()),0.0);
		} else if constexpr(dim == 3){
		glVertex3f(trans(vec.x()),trans(vec.y()),trans(vec.z()));
		} else {
			
		}
		glEnd();
	}

	inline void render_loop()
	{
		//glClearColor ( 1.0f, 1.0f, 1.0f, 1.0f );
		glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glPointSize(1);
		glLineWidth(2.5);
		
		
		glColor3f(0.0, 1.0, 0.0);
		//for(auto g: sim_setting.SimulationGrids)
		//{
		//	draw.dot(g.Position);
		//}
		
		glColor3f(0.0, 1.0, 0.0);
		for(int i = 0 ; i < sim_setting.mg.numgrids ; ++i)
		{
			//glColor3f(std::fabs(sim_setting.mg.Velocity.col(i).x()) , 0.0, std::fabs(sim_setting.mg.Velocity.col(i).y()) );
			//draw_dot<SET::dim>(sim_setting.mg.Position.col(i), 1.0 + 10.0 * std::fabs(sim_setting.mg.Velocity.col(i).y()));
		}
		
		glColor3f(1.0, 0.0, 0.0);
		for(int i = 0 ; i < sim_setting.mp.numpoints ; ++i)
		{
			glColor3f(sim_setting.mp.color[i].x(), sim_setting.mp.color[i].y(), sim_setting.mp.color[i].z());
			draw_dot<SET::dim>(sim_setting.mp.Position.col(i), 1.0);
		}
		
	}


	inline void display(void)
	{
		glClear(GL_COLOR_BUFFER_BIT);
			render_loop();
		glFlush();
	}



	inline void resize(int w, int h)
	{
		glViewport(0, 0, w, h);

		glLoadIdentity();
		glLoadIdentity();
		gluPerspective(30.0, (double)w / (double)h, 1.0, 100.0);
		gluLookAt(3.0, 4.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

		//glOrtho(-w / 200.0, w / 200.0, -h / 200.0, h / 200.0, -1.0, 1.0);
	}

	inline void idle()
	{
		integrater.integrate();
		//std::cout << "elappsed time : " <<sim_setting.iterate_num*sim_setting.dt << std::endl;
		if(sim_setting.iterate_num%io_setting.output_per_frame == 0)
		{
			glutPostRedisplay();
			std::cout << "elappsed time : " <<sim_setting.iterate_num*sim_setting.dt << std::endl;
			MPM::WRITE::show();
			//MPM::WRITE::save_dat_file();
			
			++io_setting.output_frame;
		}
		
		
		
		if(sim_setting.simulation_stop == true)
		{
			timer_main.Get_end();
			timer_main.Show_duration();
			exit(0);
		}
		
	}
}

#endif /* opengl_view_hpp */
