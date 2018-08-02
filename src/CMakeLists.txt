cmake_minimum_required(VERSION 3.2)
project(EvolutionaryLib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(EvolutionaryLib ${SOURCE_FILES} movpeaks.cpp Individual.cpp Individual.h Population.cpp Population.h Util.cpp Util.h Test.cpp Test.h Algorithm.cpp Algorithm.h Statistics.cpp Statistics.h ScoutPopulation.cpp ScoutPopulation.h ParentPopulation.cpp ParentPopulation.h)
target_link_libraries(EvolutionaryLib config++)