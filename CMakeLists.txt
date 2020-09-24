cmake_minimum_required(VERSION 3.18)
project(custom_view_select)

set(CMAKE_CXX_STANDARD 17)

include_directories(inc)

add_executable(${PROJECT_NAME} src/main.cc src/ColmapSparseInfo.cc inc/ColmapSparseInfo.hpp src/MVSNetParser.cc inc/MVSNetParser.hpp)

find_package(Eigen3 REQUIRED)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS filesystem program_options log)
target_link_libraries(${PROJECT_NAME} Boost::filesystem Boost::program_options Boost::log)

find_package( OpenCV REQUIRED )
target_include_directories( ${PROJECT_NAME} SYSTEM PUBLIC ${OpenCV_INCLUDE_DIRS} )
#target_link_libraries( ${PROJECT_NAME} ${OpenCV_LIBS})
target_link_libraries( ${PROJECT_NAME} opencv_core opencv_imgcodecs)

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
