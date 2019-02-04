cmake_minimum_required(VERSION 3.5)

project(admm)
enable_language(Fortran)

#------------------------------------------------------------
# Compiler flags
#------------------------------------------------------------

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(CMAKE_POSITION_INDEPENDENT_CODE ON)	# -fPIC
	#set(dialect "-ffree-form -std=f90 -fimplicit-none")
	set(dialect "-ffree-form -fimplicit-none")
	set(bounds "-fbounds-check")
	set(checkall "-fcheck=all")		# -Wall
	#	set(norange "-fno-range-check")	# for mkl
	set(trap "-ffpe-trap=zero,invalid,overflow,underflow")
	set(profile "-pg")
	set(O3 "-O3")
	set(O0 "-O0")
	set(Og "-Og")
	set(memflags "-g -O0")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
	#	set(dialect "-stand f90 -free -implicitnone")
	set(dialect "-free -implicitnone")
	set(bounds "-check bounds")
	# TODO: intel flag
	#	set(norange "-fno-range-check")	# for mkl
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
	# -Mallocatable=03: enable Fortran 2003 assignment semantics
	set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
	set(bounds "-C")
	# TODO: PGI flag. Not needed if not using mkl
	# set(norange "-fno-range-check")	# for mkl
	#--------------------------------------------------------
	# OpenACC
	#--------------------------------------------------------
	set(acc "-acc")
	set(gpu "-ta=tesla:cc35")
	set(cpu "-ta=multicore")
	set(omp "-mp")
	set(profile "-g -Minfo=ccff")
	set(memcheck "-g -O0")
	# Provide messages on parallel code generated
	set(feedback "-Minfo")
	# Optimizations
	set(fast "-fast")
	set(O4 "-O4")

endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
	set(dialect "-ffree")
	set(bounds "-b")
endif()

#------------------------------------------------------------
# Flags for different cmake types
#------------------------------------------------------------
set(CMAKE_Fortran_FLAGS "${dialect} ${norange} ${checkall} ${Og}")
set(CMAKE_Fortran_FLAGS_PROFILE "${dialect} ${norange} ${checkall} ${profile}")
set(CMAKE_Fortran_FLAGS_DEBUG "${dialect} ${norange} ${checkall} ${bounds} ${trap} ${profile}")
set(CMAKE_Fortran_FLAGS_MEMCHECK "${dialect} ${norange} ${checkall} ${memflags}")
set(CMAKE_Fortran_FLAGS_PERF "${dialect} ${norange} ${checkall} ${O3}")

set(CMAKE_Fortran_FLAGS_GPU "${dialect} ${acc} ${gpu} ${fast} ${feedback}")
set(CMAKE_Fortran_FLAGS_CPU "${dialect} ${acc} ${cpu} ${fast} ${feedback}")
set(CMAKE_Fortran_FLAGS_OMP "${dialect} ${acc} ${omp}")
set(CMAKE_Fortran_FLAGS_GPU_PROF "${dialect} ${acc} ${gpu} ${profile} ${feedback}")
set(CMAKE_Fortran_FLAGS_CPU_PROF "${dialect} ${acc} ${cpu} ${profile}")
set(CMAKE_Fortran_FLAGS_OMP_PROF "${dialect} ${acc} ${omp} ${profile} ${feedback}")

#------------------------------------------------------------
# Compile
#------------------------------------------------------------
file(GLOB_RECURSE sources src/*.f90)
file(GLOB_RECURSE tests test/*.f90)
file(GLOB_RECURSE external external/*.f90)
#file(GLOB_RECURSE mkl /opt/intel/compilers_and_libraries/linux/mkl/include/*.f90)

#------------------------------------------------------------
# Make static library to ADMM engine
#------------------------------------------------------------
#add_library(admm ${sources} ${external} ${mkl})
#add_library(admm ${sources} ${external})
add_library(admm STATIC ${sources} ${external})
#add_library(admm SHARED ${sources} ${external})

#------------------------------------------------------------
# USE coarrays
#------------------------------------------------------------
# Add coarray support
# if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
# 	set(coflags -coarray=shared)
# 	set(OpenCoarrays_LIBRARIES -coarray=shared)
# elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
# 	list(APPEND OpenCoarrays_LIBRARIES caf_mpi)
# 	set(coflags -fcoarray=lib)
# endif()
#------------------------------------------------------------

#------------------------------------------------------------
# UMFPACK library
#------------------------------------------------------------
# Add UMFPACK wrapper and link library
# list(APPEND UMFPACK_LIB umfpack)  
# Must have already made and installed SuiteSparse
# Can install suitesparse (archlinux)
list(APPEND UMFPACK_LIB umfpack amd cholmod suitesparseconfig)

#------------------------------------------------------------
# MKL
#------------------------------------------------------------
# MKL
#list(APPEND MKL_LIB mkl_gf_ilp64 mkl_sequential mkl_core) # mkl_gnu_thread)

#------------------------------------------------------------
# LAPACK
#------------------------------------------------------------
list(APPEND LAPACK_LIB lapack)

#------------------------------------------------------------
# HSL: MA57
#------------------------------------------------------------
list(APPEND MA57_LIB ma57 metis-4)

#------------------------------------------------------------
# Make executable for all tests
#------------------------------------------------------------
foreach(test_src ${tests})
	# cut off .o
	get_filename_component(testname ${test_src} NAME_WE)
	add_executable(${testname} ${test_src})
	# Compile options

	#--------------------------------------------------------
	# USE coarrays
	# target_compile_options(${testname} PRIVATE ${coflags})
	#--------------------------------------------------------

	# Link library to each exec
	target_link_libraries(${testname} admm)

	#--------------------------------------------------------
	# Link coarrays libs
	# target_link_libraries(${testname} ${OpenCoarrays_LIBRARIES})
	#--------------------------------------------------------
	# Link UMFPACK C original libs
	target_link_libraries(${testname} ${UMFPACK_LIB})
	# MKL
	# target_link_libraries(${testname} ${MKL_LIB})
	# Link MA57
	target_link_libraries(${testname} ${MA57_LIB})

endforeach(test_src ${tests})
