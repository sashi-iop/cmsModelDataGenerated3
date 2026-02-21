#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(cmsmodel)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# See the documentation for a guide on how to enable/disable specific components
#
find_package(Geant4 REQUIRED ui_all vis_all)

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, use our local headers, and link it to the Geant4 libraries
#
add_executable(cmsmodel serc19_cmsmodel.cc ${sources} ${headers})
target_include_directories(cmsmodel PRIVATE include)
target_link_libraries(cmsmodel ${Geant4_LIBRARIES})
target_link_libraries(cmsmodel ${ROOT_LIBRARIES})
target_compile_definitions(cmsmodel PRIVATE G4VIS_USE)
target_compile_definitions(cmsmodel PRIVATE G4UI_USE)

#----------------------------------------------------------------------------
# Build the digitisation macro as a standalone executable
#
add_executable(ecal_hcal_digitisation skeleton_ecal_hcal_digitisation.C)
target_include_directories(ecal_hcal_digitisation PRIVATE include)
target_link_libraries(ecal_hcal_digitisation ${ROOT_LIBRARIES})
target_link_libraries(ecal_hcal_digitisation ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Build the track finder as a standalone executable
#
add_executable(track_finder skelaton_track_fitting.C)
target_include_directories(track_finder PRIVATE include)
target_link_libraries(track_finder ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build cmsmodel. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(cmsmodel_SCRIPTS
  gui.mac
  icons.mac
  init_vis.mac
  run.mac
  vis.mac
)

foreach(_script ${cmsmodel_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
