cmake_minimum_required(VERSION 3.9)
project(custom_view_select)

set(CMAKE_CXX_STANDARD 17)

include_directories(inc)

add_executable(${PROJECT_NAME} src/main.cc
  src/ColmapSparseInfo.cc inc/ColmapSparseInfo.hpp
  src/MVSNetParser.cc inc/MVSNetParser.hpp
  src/Options.cc inc/Options.hpp
  src/MVSNetViewSelect.cc inc/MVSNetViewSelect.hpp
  src/ColmapViewSelect.cc inc/ColmapViewSelect.hpp)

find_package(Eigen3 REQUIRED)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS filesystem program_options log)
target_link_libraries(${PROJECT_NAME} Boost::filesystem Boost::program_options Boost::log)
target_compile_definitions(${PROJECT_NAME} PUBLIC -DBOOST_LOG_DYN_LINK)

find_package( OpenCV REQUIRED )
target_include_directories( ${PROJECT_NAME} SYSTEM PUBLIC ${OpenCV_INCLUDE_DIRS} )
#target_link_libraries( ${PROJECT_NAME} ${OpenCV_LIBS})
target_link_libraries( ${PROJECT_NAME} opencv_core opencv_imgcodecs)

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
