cmake_minimum_required(VERSION 2.8.3)
project(cans)
## Project name can now be referred to with ${PROJECT_NAME}

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  pcl_conversions
  pcl_ros
  roscpp
  sensor_msgs
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

# find_package(ImageMagick COMPONENTS Magick++ REQUIRED)

find_package(NURBS++ REQUIRED)

# find_package(GTest REQUIRED)

## "wet" packages... -> I guess this means catkin packages...
## System dependencies are found with CMake's conventions
# find_package(Boost REQUIRED COMPONENTS system)


#####################################
## Things I have tried 
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/include/nurbs++")

#set(NURBS++_INCLUDE_DIRS "/usr/local/include/nurbs++/")
#set(NURBS++_LIBRARIES "/usr/local/lib/")

#FIND_LIBRARY(NURBS_LIBRARY nurbs /usr/local/lib/)

## Add the PkgConfig package 
#find_package(PkgConfig REQUIRED)
## Get the NURBS++ module
#pkg_check_modules(NURBS++ REQUIRED nurbs++-3.0.12)
## NURBS++ is from the nurbs++-3.0.12.pc file in /usr/local/lib/pkgconfig/

#set(NURBS++_LIBRARIES ${NURBS++_LIBDIR} )
#set(NURBS++_INCLUDE_DIRS ${NURBS++_INCLUDEDIR} )


#####################################


###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if your package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
 INCLUDE_DIRS include
#  LIBRARIES nurbs_test
#  CATKIN_DEPENDS pcl_conversions pcl_ros roscpp sensor_msgs
  DEPENDS NURBS++ 
  # ImageMagick
)
## Added DEPENDS NURBS++
# name is from the nurbs++-3.0.12.pc file in /usr/local/lib/pkgconfig/

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${NURBS3_INCLUDE_DIRS}
  # ${GTEST_INCLUDE_DIRS}
  # ${ImageMagick_INCLUDE_DIRS}
  #/usr/local/include/nurbs++/
)
## Added the NURBS++ include directories here. This directory is defined in the nurbs++-3.0.12.pc file
## The _INCLUDE_DIRS is a standard suffix to denote the include path

## Declare a C++ library
## The add_library() CMake function is used to specify libraries to build. By default catkin builds shared libraries.
# add_library(${PROJECT_NAME}
#   src/mapping3D.cpp
# )

## Add cmake target dependencies of the library
## as an example, code may need to be generated before libraries
## either from message generation or dynamic reconfigure
# add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Declare a C++ executable
## With catkin_make all packages are built within a single CMake context
## The recommended prefix ensures that target names across packages don't collide
add_executable(test_mapping examples/test_mapping.cpp src/mapping3D.cpp src/object3D.cpp src/splitSurface.cpp)
add_executable(crv_test examples/curveTesting.cpp src/mapping3D.cpp src/object3D.cpp src/splitSurface.cpp)
add_executable(srf_test examples/surfaceJoining.cpp src/mapping3D.cpp src/object3D.cpp src/splitSurface.cpp)
add_executable(pcl_test examples/pointCloudTesting.cpp src/mapping3D.cpp src/object3D.cpp src/splitSurface.cpp)
add_executable(localize_test examples/localizationTesting.cpp)
add_executable(split_test examples/testingSplitSurface.cpp src/splitSurface.cpp)


# add_definitions( -DMAGICKCORE_QUANTUM_DEPTH=16 )
# add_definitions( -DMAGICKCORE_HDRI_ENABLE=0 )

## Rename C++ executable without prefix
## The above recommended prefix causes long target names, the following renames the
## target back to the shorter version for ease of user use
## e.g. "rosrun someones_pkg node" instead of "rosrun someones_pkg someones_pkg_node"
# set_target_properties(${PROJECT_NAME}_node PROPERTIES OUTPUT_NAME node PREFIX "")

## Add cmake target dependencies of the executable
## same as for the library above
# add_dependencies(${PROJECT_NAME}_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Specify libraries to link a library or executable target against
target_link_libraries(test_mapping ${catkin_LIBRARIES} ${NURBS3_LIBRARIES}) 
target_link_libraries(crv_test ${catkin_LIBRARIES} ${NURBS3_LIBRARIES}) 
target_link_libraries(srf_test ${catkin_LIBRARIES} ${NURBS3_LIBRARIES}) 
target_link_libraries(pcl_test ${catkin_LIBRARIES} ${NURBS3_LIBRARIES}) 
target_link_libraries(localize_test ${catkin_LIBRARIES} ${NURBS3_LIBRARIES}) 
target_link_libraries(split_test ${catkin_LIBRARIES} ${NURBS3_LIBRARIES}) 

#############
## Install ##
#############
# Probably don't care too much about this yet - don't want to install quite yet

# all install targets should use catkin DESTINATION variables
# See http://ros.org/doc/api/catkin/html/adv_user_guide/variables.html

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
# install(PROGRAMS
#   scripts/my_python_script
#   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark executables and/or libraries for installation
# install(TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_node
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark cpp header files for installation
# install(DIRECTORY include/${PROJECT_NAME}/
#   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
#   FILES_MATCHING PATTERN "*.h"
#   PATTERN ".svn" EXCLUDE
# )

## Mark other files for installation (e.g. launch and bag files, etc.)
# install(FILES
#   # myfile1
#   # myfile2
#   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )

#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_nurbs.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()
if (CATKIN_ENABLE_TESTING)
  catkin_add_gtest(test_mapping3D tests/test_mapping3D.cpp src/mapping3D.cpp src/object3D.cpp src/splitSurface.cpp)
  target_link_libraries(test_mapping3D ${catkin_LIBRARIES} ${NURBS3_LIBRARIES})
  catkin_add_gtest(test_object3D tests/test_object3D.cpp src/object3D.cpp)
  target_link_libraries(test_object3D ${catkin_LIBRARIES} ${NURBS3_LIBRARIES})
  catkin_add_gtest(test_splitSurface tests/test_splitSurface.cpp src/splitSurface.cpp)
  target_link_libraries(test_splitSurface ${catkin_LIBRARIES} ${NURBS3_LIBRARIES})
endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
