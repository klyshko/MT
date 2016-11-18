project(MT)

cmake_minimum_required(VERSION 2.8)

include_directories(.)


set(CMAKE_CXX_FLAGS -fPIC)
set(CMAKE_C_FLAGS -fPIC)

option(MORSE  "Use morse potential." ON)
option(READ_FROM_DCD  "Read coordinates from dcd" OFF)
option(WITH_MPI "Build with MPI support." OFF)
option(AVERAGE_LJ "Output average LJ force and st deviation" OFF)
option(REGULARIZATION "Regularization of LJ force. Making potential finite" OFF)


if(${WITH_MPI})
	find_package(MPI REQUIRED)
	include_directories(${MPI_INCLUDE_PATH})
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D USE_MPI")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D USE_MPI")
endif()

if(${REGULARIZATION})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D REGULARIZATION")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D REGULARIZATION")
endif()
										
if(${AVERAGE_LJ})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D AVERAGE_LJ")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D AVERAGE_LJ")
endif()


if(${MORSE})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D MORSE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D MORSE")
endif()

if(${READ_FROM_DCD})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D READ_FROM_DCD")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D READ_FROM_DCD")
endif()

# Find CUDA
if (NOT DEFINED CUDA_ARCH)
  set(CUDA_ARCH sm_35)
endif(NOT DEFINED CUDA_ARCH)

if (NOT DEFINED CUDA_SEPARABLE_COMPILATION)
  set(CUDA_SEPARABLE_COMPILATION ON)
endif(NOT DEFINED CUDA_SEPARABLE_COMPILATION)

set(CUDA_ARCH "${CUDA_ARCH}" CACHE STRING "GPU Architecture" FORCE)
find_package(CUDA QUIET REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${CUDA_SDK_ROOT_DIR}/common/inc)

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch ${CUDA_ARCH} -use_fast_math)
#set(CUDA_NVCC_FLAGS -arch ${CUDA_ARCH} --ptxas-options=-v ${CUDA_NVCC_FLAGS})

# Disable OpenMP
set(ENABLE_OMP OFF)

option (CUDA_ENFORCE_HOST_COMPILER "Force nvcc to use the same compiler used to compile .c(pp) files insted of gcc/g++" OFF)
if (${CUDA_ENFORCE_HOST_COMPILER})
  set (CMAKE_GCC_TEMP_DIR "CMakeGCC")
  file(MAKE_DIRECTORY ${CMAKE_GCC_TEMP_DIR})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_C_COMPILER} ${CMAKE_GCC_TEMP_DIR}/gcc)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CXX_COMPILER} ${CMAKE_GCC_TEMP_DIR}/g++)
  set(CUDA_NVCC_FLAGS -ccbin ${CMAKE_GCC_TEMP_DIR} ${CUDA_NVCC_FLAGS})
endif(${CUDA_ENFORCE_HOST_COMPILER})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D CUDA")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D CUDA")

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo." FORCE)
endif()


set(SOURCES
		src/dcdio.cpp
		src/xyzio.cpp
    		src/globals.cpp
		src/main.cpp
		src/pdbio.cpp
		src/configreader.cpp
		src/wrapper.cpp
		src/timer.cpp
		src/parameters.cpp
    		src/compute_cuda.cu
    		src/updater.cpp
    		src/preparator.cpp
    		src/bdhitea.cu
    		src/bdhitea_kernel.cu
    		src/HybridTaus.cu

)

cuda_add_executable(mt ${SOURCES})

target_link_libraries(mt ${CUDA_LIBRARIES} ${MPI_LIBRARIES})
