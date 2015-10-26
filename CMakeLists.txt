project(MT)

cmake_minimum_required(VERSION 2.8)

include_directories(.)


set(CMAKE_CXX_FLAGS -fPIC)
set(CMAKE_C_FLAGS -fPIC)


option(OUTPUT_EN  "Output energy" ON)
option(ANGLES  "Use complex representations of angle interactions" OFF)
option(SIMPLEWELL  "Temporary and final feature" OFF)
option(OUTPUT_FORCE "Output force" ON)
option(SIGMOID  "Use sigmoid function unstead of exponential or harmonic." OFF)
option(DEXP  "Use double-Gaussian potential." OFF)
option(EXP  "Use exponential potential." OFF)
option(MORSE  "Use morse potential." ON)
option(BARR  "Add barrier" OFF)
option(LJ  "Use LJ potential" ON)
option(READ_FROM_DCD  "Read coordinates from dcd" OFF)
option(DOUBLEWELL  "Use ax^4+bx^2 instead of ax^2 for angles." OFF)
option(REDUCE_TO_2D  "Freeze any movement along Y, psi, phi to make this model behave like 2D." OFF)
option(WITH_MPI "Build with MPI support." OFF)
option(REPULSIVE "Turn on repulsive potential." OFF)
option(ASSEMBLY "Turn on assembly mode." OFF)
option(AVERAGE_LJ "Output average LJ force and st deviation" OFF)
option(REGULARIZATION "Regularization of LJ force. Making potential finite" OFF)
option(MT_LENGTH "Write MT length to file" ON)

if(${WITH_MPI})
	find_package(MPI REQUIRED)
	include_directories(${MPI_INCLUDE_PATH})
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D USE_MPI")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D USE_MPI")
endif()

if(${MT_LENGTH})
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D MT_LENGTH")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D MT_LENGTH")
endif()

if(${REGULARIZATION})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D REGULARIZATION")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D REGULARIZATION")
endif()
										

if(${OUTPUT_EN})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D OUTPUT_EN")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D OUTPUT_EN")
endif()

if(${AVERAGE_LJ})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D AVERAGE_LJ")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D AVERAGE_LJ")
endif()

if(${BARR})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D BARR")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D BARR")
endif()
if(${EXP})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D EXP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D EXP")
endif()
if(${MORSE})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D MORSE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D MORSE")
endif()
if(${ANGLES})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D ANGLES")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D ANGLES")
endif()
if(${OUTPUT_FORCE})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D OUTPUT_FORCE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D OUTPUT_FORCE")
endif()
if(${READ_FROM_DCD})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D READ_FROM_DCD")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D READ_FROM_DCD")
endif()
if(${SIGMOID})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D SIGMOID")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D SIGMOID")
endif()
if(${DEXP})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D DEXP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D DEXP")
endif()
if(${REDUCE_TO_2D})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D REDUCE_TO_2D")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D REDUCE_TO_2D")
endif()
if(${DOUBLEWELL})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D DOUBLEWELL")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D DOUBLEWELL")
endif()
if(${SIMPLEWELL})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D SIMPLEWELL")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D SIMPLEWELL")
endif()
if(${LJ})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D LJ_on")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D LJ_on")
endif()
if(${REPULSIVE})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D REPULSIVE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D REPULSIVE")
endif()
if(${ASSEMBLY})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D ASSEMBLY")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D ASSEMBLY")
endif()

# Find CUDA
if (NOT DEFINED CUDA_ARCH)
  set(CUDA_ARCH sm_20)
endif(NOT DEFINED CUDA_ARCH)
set(CUDA_ARCH "${CUDA_ARCH}" CACHE STRING "GPU Architecture" FORCE)
find_package(CUDA QUIET REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${CUDA_SDK_ROOT_DIR}/common/inc)
set(CUDA_NVCC_FLAGS -arch ${CUDA_ARCH} --ptxas-options=-v -use_fast_math ${CUDA_NVCC_FLAGS})
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
		src/main.cpp
		src/pdbio.cpp
		src/configreader.cpp
		src/wrapper.cpp
		src/timer.cpp
		src/parameters.cpp
        src/compute_cuda.cu
)

cuda_add_executable(mt ${SOURCES})

target_link_libraries(mt ${CUDA_LIBRARIES} ${MPI_LIBRARIES})
