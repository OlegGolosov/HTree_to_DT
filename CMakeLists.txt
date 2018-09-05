cmake_minimum_required(VERSION 2.4 FATAL_ERROR)

project(HTree_TO_DT)

set(CMAKE_BUILD_TYPE debug)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH $ENV{DATATREE_HOME})

find_package(ROOT REQUIRED COMPONENTS RIO Tree)
find_package(DataTree REQUIRED)

set(CMAKE_CXX_FLAGS "-Wall")

include(${ROOT_USE_FILE})

include_directories(
    $ENV{HADDIR}/include
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
    ${DataTree_INCLUDE_DIR}
    )
		
link_directories(
    $ENV{HADDIR}/lib
    ${CMAKE_SOURCE_DIR}
    ${DataTree_LIBRARY_DIR}
    )
set (HADES_LIBS
		libAlignment.so     
		libHydra.so        
		libMdcTrackD.so  
		libOra.so       
		libPionTracker.so  
		libRpc.so         
		libStart.so
		libDst.so           
		libKalman.so       
		libMdcTrackG.so  
		libOraSim.so    
		libQA.so           
		libShower.so      
		libTof.so
		libEmc.so           
		libMdc.so          
		libMdcUtil.so   
		libOraUtil.so   
		libRevt.so         
		libShowerUtil.so  
		libTools.so
		libEventDisplay.so  
		libMdcGarfield.so  
		libOnline.so     
		libParticle.so  
		libRich.so         
		libSimulation.so  
		libWall.so
		)
add_definitions(${ROOT_CXX_FLAGS})

add_executable(HTree_to_DT src/HADES_constants.h src/mhwalldivider.h src/HTree_to_DT.cxx)
target_link_libraries(HTree_to_DT ${ROOT_LIBRARIES} ${HADES_LIBS} DataTree)
