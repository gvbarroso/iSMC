# CMake script for iSMC
# Authors:
#   Julien Dutheil
# Created: 16/04/2018

SET(ismc_CPP
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

SET(ismc_mapper_CPP
  mapper.cpp
)

SET(ismc_mapper2_CPP
  mapper2.cpp
)

ADD_EXECUTABLE (ismc ${ismc_CPP})
ADD_EXECUTABLE (ismc_mapper ${ismc_mapper_CPP})
ADD_EXECUTABLE (ismc_mapper2 ${ismc_mapper2_CPP})

SET (ismc-targets ismc ismc_mapper ismc_mapper2)

FOREACH (target ${ismc-targets})
  # Link (static or shared)
  IF (BUILD_STATIC)
    SET_TARGET_PROPERTIES (${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
    TARGET_LINK_LIBRARIES (${target} ${BPP_LIBS_STATIC} ${ZIPHMM_LIBS_STATIC})
    TARGET_LINK_LIBRARIES (${target} ${LIBS})
    TARGET_LINK_LIBRARIES (${target} Threads::Threads)
  ELSE (BUILD_STATIC)
    TARGET_LINK_LIBRARIES (${target} ${BPP_LIBS_SHARED} ${ZIPHMM_LIBS_SHARED})
    TARGET_LINK_LIBRARIES (${target} ${LIBS})
    TARGET_LINK_LIBRARIES (${target} Threads::Threads)
    SET_TARGET_PROPERTIES (${target} PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  ENDIF (BUILD_STATIC)
ENDFOREACH (target)

INSTALL(TARGETS ${ismc-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})

