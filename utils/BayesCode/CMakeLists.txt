cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(BayesCode CXX)

# Compilation options
option(COVERAGE_MODE "For coverage mode using g++ " OFF) #OFF by default
option(TINY "Compiling a subset of programs without MPI " OFF) #OFF by default
option(DEBUG_MODE "Debug mode (with asserts and such) " OFF) #OFF by default
set(ALWAYSON_CXX_FLAGS "--std=c++14 -Wfatal-errors ")
if (COVERAGE_MODE)
    set(CMAKE_CXX_FLAGS "-O0 -fprofile-arcs -ftest-coverage ${ALWAYSON_CXX_FLAGS}") # coverage mode
    message("-- INFO: Compiling in coverage mode.\n-- INFO: flags are: " ${CMAKE_CXX_FLAGS})
else ()
    if (DEBUG_MODE)
        set(CMAKE_CXX_FLAGS "-O2 -g ${ALWAYSON_CXX_FLAGS}") # debug mode
        message("-- INFO: Compiling in debug mode.\n-- INFO: flags are: " ${CMAKE_CXX_FLAGS})
    else ()
        set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG ${ALWAYSON_CXX_FLAGS}") # release mode
        message("-- INFO: Compiling in release mode.\n-- INFO: flags are: " ${CMAKE_CXX_FLAGS})
    endif (DEBUG_MODE)
endif (COVERAGE_MODE)
unset(COVERAGE_MODE CACHE)
unset(DEBUG_MODE CACHE)

if (TINY)
    message("-- INFO: Compiling without MPI.")
else ()
    message("-- INFO: Compiling with MPI.")
    # MPI
    find_package(MPI REQUIRED)
    include_directories(${MPI_INCLUDE_PATH})
endif (TINY)

include_directories("src")
include_directories("src/lib")
include_directories("utils")

# Static libs
add_library(tree_lib STATIC "src/tree/implem.cpp")

set(BAYESCODE_LIB
        src/lib/AAMutSelOmegaCodonSubMatrix.cpp
        src/lib/AAMutSelNeCodonMatrixBidimArray.cpp
        src/lib/BranchSitePath.cpp
        src/lib/Chain.cpp
        src/lib/Chrono.cpp
        src/lib/CodonSequenceAlignment.cpp
        src/lib/CodonStateSpace.cpp
        src/lib/CodonSubMatrix.cpp
        src/lib/GTRSubMatrix.cpp
        src/lib/PhyloProcess.cpp
        src/lib/Random.cpp
        src/lib/RandomGeneSample.cpp
        src/lib/Sample.cpp
        src/lib/SequenceAlignment.cpp
        src/lib/StateSpace.cpp
        src/lib/SubMatrix.cpp
        src/lib/TaxonSet.cpp
        src/lib/Chronogram.cpp
        src/lib/MultivariateProcess.cpp
        src/lib/ScaledMutationRateCompound.cpp
        src/lib/ScatterSuffStat.cpp
        src/lib/PolyData.cpp
        src/lib/PolyProcess.cpp
        src/lib/PolySuffStat.cpp
        src/lib/PoissonRandomField.cpp
        )
add_library(bayescode_lib STATIC ${BAYESCODE_LIB})

set(BASE_LIBS
        bayescode_lib
        tree_lib
        )

# Mutation-Selection (mutselomega)
add_executable(mutselomega "src/AAMutSelMultipleOmega.cpp")
target_link_libraries(mutselomega ${BASE_LIBS})

# Read Mutation-Selection (readmutselomega)
add_executable(readmutselomega "src/ReadAAMutSelMultipleOmega.cpp")
target_link_libraries(readmutselomega ${BASE_LIBS})

# Node Ne mutsel (nodemutsel)
add_executable(nodemutsel "src/DatedNodeMutSel.cpp")
target_link_libraries(nodemutsel ${BASE_LIBS})

# Read Node Ne mutsel (readnodemutsel)
add_executable(readnodemutsel "src/ReadDatedNodeMutSel.cpp")
target_link_libraries(readnodemutsel ${BASE_LIBS})

# Read Node Brownian (readnodetraitsmutsel)
add_executable(readnodetraits "src/ReadDatedNode.cpp")
target_link_libraries(readnodetraits ${BASE_LIBS})

# Node Brownian (nodetraitsmutsel)
add_executable(nodetraits "src/DatedNode.cpp")
target_link_libraries(nodetraits ${BASE_LIBS})

# Node omega (nodeomega)
add_executable(nodeomega "src/DatedNodeOmega.cpp")
target_link_libraries(nodeomega ${BASE_LIBS})

# Read Node omega (readnodeomega)
add_executable(readnodeomega "src/ReadDatedNodeOmega.cpp")
target_link_libraries(readnodeomega ${BASE_LIBS})

if (TINY)
    message("-- INFO: Compiling uniquely mutselomega, nodetraits, nodemutsel and nodeomega.")
else ()
    # Single omega (globom)
    add_executable(globom "src/SingleOmega.cpp")
    target_link_libraries(globom ${BASE_LIBS})

    # Read single omega (readglobom)
    add_executable(readglobom "src/ReadSingleOmega.cpp")
    target_link_libraries(readglobom ${BASE_LIBS})

    # Codon M2a
    add_executable(codonm2a "src/CodonM2a.cpp" "src/CodonM2aModel.cpp")
    target_link_libraries(codonm2a ${BASE_LIBS})

    # Read Codon M2a
    add_executable(readcodonm2a "src/ReadCodonM2a.cpp" "src/CodonM2aModel.cpp")
    target_link_libraries(readcodonm2a ${BASE_LIBS})

    # Mutation-Selection (aamutsel)
    add_executable(aamutsel "src/AAMutSelDSBDPOmega.cpp")
    target_link_libraries(aamutsel ${BASE_LIBS})

    # Read Mutation-Selection (readaamutsel)
    add_executable(readaamutsel "src/ReadAAMutSelDSBDPOmega.cpp")
    target_link_libraries(readaamutsel ${BASE_LIBS})

    # Mutation-Selection (mutseldm5)
    add_executable(mutseldm5 "src/AAMutSelDM5.cpp")
    target_link_libraries(mutseldm5 ${BASE_LIBS})

    # Read Mutation-Selection (readmutseldm5)
    add_executable(readmutseldm5 "src/ReadAAMutSelDM5.cpp")
    target_link_libraries(readmutseldm5 ${BASE_LIBS})

    # Diffsel double sparse (diffseldsparse)
    add_executable(diffseldsparse "src/DiffSelDoublySparse.cpp")
    target_link_libraries(diffseldsparse ${BASE_LIBS})

    # Read diffsel double sparse (readdiffseldsparse)
    # add_executable(readdiffseldsparse "src/ReadDiffSelDoublySparse.cpp")
    # target_link_libraries(readdiffseldsparse ${BASE_LIBS})

    # # Multi-gene single omega (multigeneglobom)
    # add_executable(multigeneglobom "src/MultiGeneSingleOmega.cpp")
    # target_link_libraries(multigeneglobom ${BASE_LIBS} ${MPI_LIBRARIES})

    # # Multi Gene Codon M2a
    # add_executable(multigenecodonm2a "src/MultiGeneCodonM2a.cpp" "src/CodonM2aModel.cpp")
    # target_link_libraries(multigenecodonm2a ${BASE_LIBS} ${MPI_LIBRARIES})

    # tests
    add_executable(tree_test "src/tree/test.cpp")
    target_link_libraries(tree_test tree_lib)
    # MPI
    add_executable(mpi_par_test "src/mpi_components/test/mpi_test.cpp")
    target_link_libraries(mpi_par_test ${MPI_LIBRARIES})

    add_executable(all_tests "src/all_tests.cpp")
    target_link_libraries(all_tests ${MPI_LIBRARIES})
endif (TINY)
unset(TINY CACHE)