cmake_minimum_required(VERSION 3.5)

project(Boundary)

find_package(Eigen3 REQUIRED)

add_executable(${PROJECT_NAME} boundary_angle.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE
    Eigen3::Eigen)