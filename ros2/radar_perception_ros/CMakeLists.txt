cmake_minimum_required(VERSION 3.20)
project(radar_perception_ros LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(visualization_msgs REQUIRED)
find_package(rosgraph_msgs REQUIRED)
find_package(tf2_ros REQUIRED)
find_package(radar_perception_msgs REQUIRED)
find_package(radar_perception 0.1 REQUIRED)  # core library: radar::core

set(ROS_DEPS rclcpp std_msgs geometry_msgs sensor_msgs nav_msgs visualization_msgs
             rosgraph_msgs tf2_ros radar_perception_msgs)

foreach(node sim_source_node radar_perception_node ego_estimator_node)
  add_executable(${node} src/${node}.cpp)
  target_compile_options(${node} PRIVATE -Wall -Wextra -Wpedantic)
  target_link_libraries(${node} radar::core)
  ament_target_dependencies(${node} ${ROS_DEPS})
  install(TARGETS ${node} DESTINATION lib/${PROJECT_NAME})
endforeach()

install(DIRECTORY launch config DESTINATION share/${PROJECT_NAME})

ament_package()
