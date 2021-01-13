//
//  timer.h
//  MPM2D
//
//  Created by Kentaro Nagasawa on 2020/12/13.
//  Copyright © 2020年 Kentaro Nagasawa. All rights reserved.
//

#ifndef timer_h
#define timer_h

#include <time.h>

class timer{
	clock_t start;
	clock_t end;
	public:
		void Get_start();
		void Get_end();
		void Show_duration();
};

inline void timer::Get_start(){
	start = clock();
	std::cout << "set start clock: " << start << std::endl;
}
inline void timer::Get_end(){
	end = clock();
	std::cout << "set end clock: " << start << std::endl;
}
inline void timer::Show_duration(){
	std::cout << "durations = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
}

#endif /* timer_h */
