cmake_minimum_required(VERSION 3.12)


project(LxSpatialOps VERSION 0.1.0)

#set(FORCE_STATIC 1)
if(DEFINED FORCE_STATIC )
	message(STATUS "NOTICE: **** Static build for templateMatchingAlignment *****")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
	add_definitions(-DIO_DATA_STATIC -DLXGEO_FACTORY_STATIC)
else()
	message(STATUS "NOTICE: **** Dynamic build templateMatchingAlignment *****")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 20)

add_subdirectory("$ENV{LX_GEO_REPOS_ROOT}/LxGeoCppLibs/libs/shared" LxGeoFactoryShared)
add_subdirectory(lib)