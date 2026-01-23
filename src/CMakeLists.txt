# CMake script for iSMC
# Authors:
#   Julien Dutheil
# Created: 16/04/2018

set(ismc_CPP
    AbstractParametrizedMatrix.cpp
    BackupListenerOv.cpp
    BaumWelch.cpp
    HmmStatesLibrary.cpp
    MarkovModulatedSmc.cpp
    MatrixOptimizationFunction.cpp
    MmPsmc.cpp
    MmSmcEmissionProbabilities.cpp
    MmSmcTransitionProbabilities.cpp
    MultipleMmPsmc.cpp
    ParameterCategoryTransitions.cpp
    PolymorphismData.cpp
    Psmc.cpp
    RelativeStopCondition.cpp
    SequentiallyMarkovCoalescent.cpp
    SmcDecodingWrapper.cpp
    SmcEmissionProbabilities.cpp
    SmcOptimizationWrapper.cpp
    SmcTransitionProbabilities.cpp
    Splines.cpp
    SplinesModel.cpp
    Vcf.cpp
    Global.cpp
    ismc.cpp
)

set(ismc_mapper_CPP mapper.cpp)

set(ismc_mapper2_CPP mapper2.cpp)

add_executable(ismc ${ismc_CPP})
add_executable(ismc_mapper ${ismc_mapper_CPP})
add_executable(ismc_mapper2 ${ismc_mapper2_CPP})

set(ismc-targets ismc ismc_mapper ismc_mapper2)

foreach(target ${ismc-targets})
    # Link (static or shared)
    if(BUILD_STATIC)
        target_link_libraries(${target} ${BPP_LIBS_STATIC} ${ZIPHMM_LIBS_STATIC} ${LIBS} Threads::Threads)
        set_target_properties(${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
	target_link_options(${target} PRIVATE -static-libgcc -static-libstdc++)
    else(BUILD_STATIC)
	target_link_libraries(${target} ${BPP_LIBS_SHARED} ${ZIPHMM_LIBS_SHARED} ${LIBS} Threads::Threads)
        set_target_properties(
            ${target}
            PROPERTIES
                INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
        )
    endif(BUILD_STATIC)
endforeach(target)

install(TARGETS ${ismc-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
