# 3.2.1 to allow set_property( TARGET target PROPERTY CXX_STANDARD 11 ) and similar
cmake_minimum_required( VERSION 3.2.1 )

option( USE_OPENGL "View OpneGL" ON )

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE debug)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")


add_executable(Main2D
  ./src/2DMain/main.cpp
)

add_executable(Main3D
  ./src/3DMain/main.cpp
)

if( USE_OPENGL )

  target_compile_options(Main2D
    PRIVATE
      -DUSEGL
    )
    
  target_compile_options(Main3D
    PRIVATE
      -DUSEGL
    )
  set( CMAKE_CXX_FLAGS "-framework GLUT -framework OpenGL")

endif()


target_compile_options(Main2D
  PRIVATE
    -Wall
    -O3
  )

target_include_directories(Main2D
  PRIVATE
    /usr/local/include
    /usr/local/Cellar/libomp/11.0.0/include
    /Users/kn/Downloads/boost_1_75_0
    /Users/kn/Documents/PhySimCode/xproj/MPM2D/include
)
target_link_directories(Main2D
  PRIVATE
    /usr/local/Cellar/libomp/11.0.0/lib
)

target_compile_options(Main3D
  PRIVATE
    -Wall
    -O3
  )

target_include_directories(Main3D
  PRIVATE
    /usr/local/include
    /usr/local/Cellar/libomp/11.0.0/include
    /Users/kn/Downloads/boost_1_75_0
    /Users/kn/Documents/PhySimCode/xproj/MPM2D/include
)

target_link_directories(Main3D
  PRIVATE
    /usr/local/Cellar/libomp/11.0.0/lib
)

