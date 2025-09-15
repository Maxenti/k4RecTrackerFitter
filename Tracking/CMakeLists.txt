# Tracking/CMakeLists.txt
set(PackageName Tracking)
project(${PackageName})

find_package(Gaudi REQUIRED)       # Gaudi::GaudiKernel
# Try to find a GaudiAlg target; make it optional
set(GAUDI_ALG_TGT "")
if (TARGET Gaudi::GaudiAlg)
  set(GAUDI_ALG_TGT Gaudi::GaudiAlg)
elseif (TARGET Gaudi::GaudiAlgLib)
  set(GAUDI_ALG_TGT Gaudi::GaudiAlgLib)
elseif (TARGET GaudiAlgLib)
  set(GAUDI_ALG_TGT GaudiAlgLib)
endif()
message(STATUS "GaudiAlg target used: '${GAUDI_ALG_TGT}' (empty means none)")

include(ExternalData)
list(APPEND ExternalData_URL_TEMPLATES
  "https://key4hep.web.cern.ch:443/testFiles/k4RecTracker/%(hash)"
)

find_package(podio REQUIRED COMPONENTS podioRootIO)
find_package(ROOT  REQUIRED COMPONENTS RIO Tree)

find_package(MarlinUtil REQUIRED)
list(APPEND CMAKE_MODULE_PATH $ENV{PANDORAPFA}/cmakemodules)
foreach(pkg MarlinUtil PandoraSDK)
  if (${pkg}_FOUND)
    include_directories(${${pkg}_INCLUDE_DIRS})
    link_libraries(${${pkg}_LIBRARIES})
    add_definitions(${${pkg}_DEFINITIONS})
  endif()
endforeach()

# Collect sources (but DO NOT add GenFit2 files unconditionally)
file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cpp
  ${PROJECT_SOURCE_DIR}/components/*.cpp
)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

# --- ONNX Runtime + Torch are already used in your code
# (Assuming imported targets exist in your env; otherwise adjust)
gaudi_add_module(${PackageName}
  SOURCES ${sources}
  LINK
    Gaudi::GaudiKernel
    ${GAUDI_ALG_TGT}
    EDM4HEP::edm4hep
    k4FWCore::k4FWCore
    k4FWCore::k4Interface
    extensionDict
    onnxruntime::onnxruntime
    torch
    DD4hep::DDRec
    DD4hep::DDCore
    ${DELPHES_LIBRARIES}
)

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${MarlinUtil_INCLUDE_DIRS}
  ${DELPHES_INCLUDE_DIRS}
)

set_target_properties(${PackageName} PROPERTIES PUBLIC_HEADER "${headers}")

# Copy test scripts into the build tree for convenience
file(GLOB scripts ${PROJECT_SOURCE_DIR}/test/*.py)
file(COPY ${scripts} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test)

# ---------- OPTIONAL: GenFit2 ----------
option(USE_GENFIT2 "Build GenFit2-based fitter(s)" OFF)

if(USE_GENFIT2)
  # Expect the user to set GENFIT_DIR to the GenFit2 install prefix
  #   GENFIT_DIR/include/genfit/RKTrackRep.h
  #   GENFIT_DIR/lib(64)/libgenfit2.so
  set(_GF_HINTS
    $ENV{GENFIT_DIR}
    ${GENFIT_DIR}
  )
  find_path(GENFIT2_INCLUDE_DIR genfit/RKTrackRep.h
    HINTS ${_GF_HINTS}
    PATH_SUFFIXES include
  )
  find_library(GENFIT2_LIBRARY NAMES genfit2
    HINTS ${_GF_HINTS}
    PATH_SUFFIXES lib64 lib
  )

  if(NOT GENFIT2_INCLUDE_DIR OR NOT GENFIT2_LIBRARY)
    message(FATAL_ERROR "USE_GENFIT2=ON but GenFit2 not found. Set GENFIT_DIR to your GenFit2 install (with include/ and lib*/).")
  endif()

  message(STATUS "GenFit2 include: ${GENFIT2_INCLUDE_DIR}")
  message(STATUS "GenFit2 library: ${GENFIT2_LIBRARY}")

  # Add the component only when GenFit2 is available
  target_sources(${PackageName} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/components/GenFit2DCHFitter.cpp
  )
  target_include_directories(${PackageName} PRIVATE ${GENFIT2_INCLUDE_DIR})
  target_link_libraries(${PackageName} PRIVATE ${GENFIT2_LIBRARY})
  target_compile_definitions(${PackageName} PRIVATE HAVE_GENFIT2=1)

  # Optional standalone helper/tool that also uses GenFit2
  add_executable(fit_ggtf_tracks ${CMAKE_SOURCE_DIR}/Tracking/fitter/fit_ggtf_tracks.cpp)
  target_include_directories(fit_ggtf_tracks PRIVATE ${GENFIT2_INCLUDE_DIR})
  target_link_libraries(fit_ggtf_tracks
    PRIVATE
      ${GENFIT2_LIBRARY}
      EDM4HEP::edm4hep
      podio::podioRootIO
      ROOT::Tree ROOT::RIO
      k4FWCore::k4FWCore
      k4FWCore::k4Interface
      extensionDict
  )
  target_compile_features(fit_ggtf_tracks PRIVATE cxx_std_17)
else()
  message(STATUS "GenFit2 disabled (USE_GENFIT2=OFF). GenFit2DCHFitter.cpp is not compiled.")
endif()

# ---------- Install ----------
install(TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}" COMPONENT dev
)

# Install test scripts, as before
install(FILES ${scripts} DESTINATION test)

# ---------- Existing tests ----------
SET(test_name "test_TracksFromGenParticles")
ADD_TEST(NAME ${test_name} COMMAND sh +x test_TracksFromGenParticles.sh)
set_test_env(${test_name})
set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/Tracking/test/testTracksFromGenParticles")

SET(test_name "test_TrackFinder")
ExternalData_Add_Test(test_TrackFinder
  NAME ${test_name}
  COMMAND sh +x Tracking/test/testTrackFinder/test_trackFinder.sh DATA{${CMAKE_SOURCE_DIR}/Tracking/test/inputFiles/SimpleGatrIDEAv3o1.onnx})
set_test_env(${test_name})
set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
ExternalData_Add_Target(${test_name})

SET(test_name "Noddsim_track_test")
ExternalData_Add_Test(Noddsim_track_test
  NAME ${test_name}
  COMMAND sh +x Tracking/test/testTrackFinder/Noddsim_track_test.sh DATA{${CMAKE_SOURCE_DIR}/Tracking/test/inputFiles/SimpleGatrIDEAv3o1.onnx})
set_test_env(${test_name})
set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
ExternalData_Add_Target(${test_name})
