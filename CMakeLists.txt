
set(DIL2AL_SOURCES

    ${SRC_DIR}/dil2al/regex-gnu.cc
    ${SRC_DIR}/dil2al/regex-gnu.h
    ${SRC_DIR}/dil2al/BigRegex.cc
    ${SRC_DIR}/dil2al/BigRegex.hh
    ${SRC_DIR}/dil2al/BigString.cc
    ${SRC_DIR}/dil2al/BigString.hh
    ${SRC_DIR}/dil2al/rx.h


)

set(GEOMETRY_SOURCES

    ${SRC_DIR}/geometry/point.c
    ${SRC_DIR}/geometry/point.h
    ${SRC_DIR}/geometry/Point2D.cc
    ${SRC_DIR}/geometry/Point2D.hh
    ${SRC_DIR}/geometry/Point3D.cc
    ${SRC_DIR}/geometry/Point3D.hh
    ${SRC_DIR}/geometry/ptest.c
    ${SRC_DIR}/geometry/vector.c
    ${SRC_DIR}/geometry/vector.h
    ${SRC_DIR}/geometry/Vector2D.cc
    ${SRC_DIR}/geometry/Vector2D.hh
    ${SRC_DIR}/geometry/Vector3D.cc
    ${SRC_DIR}/geometry/Vector3D.hh
    ${SRC_DIR}/geometry/vtest.c

)


# Configure Sources
set(MAIN_SOURCES 

    ${SRC_DIR}/common_entry.cc

    ${SRC_DIR}/axon_direction_model.cc
    ${SRC_DIR}/axon_direction_model.hh

    ${SRC_DIR}/BigRegex.cc
    ${SRC_DIR}/BigRegex.hh
    ${SRC_DIR}/BigString.cc
    ${SRC_DIR}/BigString.hh

    ${SRC_DIR}/branching_models.cc
    ${SRC_DIR}/branching_models.hh

    ${SRC_DIR}/Catacomb_Object.cc
    ${SRC_DIR}/Catacomb_Object.hh

    ${SRC_DIR}/Color_Table.cc
    ${SRC_DIR}/Color_Table.hh
    ${SRC_DIR}/Command_Line_Parameters.cc
    ${SRC_DIR}/Command_Line_Parameters.hh
    ${SRC_DIR}/Command.cc
    ${SRC_DIR}/Command.hh
    ${SRC_DIR}/common.h
    ${SRC_DIR}/Connection_Statistics.cc
    ${SRC_DIR}/Connection_Statistics.hh
    ${SRC_DIR}/connection.cc
    ${SRC_DIR}/connection.hh
    ${SRC_DIR}/dendritic_growth_model.cc
    ${SRC_DIR}/dendritic_growth_model.hh
    ${SRC_DIR}/diagnostic.cc
    ${SRC_DIR}/diagnostic.hh
    ${SRC_DIR}/environment_physics.cc
    ${SRC_DIR}/environment_physics.hh
    ${SRC_DIR}/event.cc
    ${SRC_DIR}/event.hh
    ${SRC_DIR}/fibre_structure.cc
    ${SRC_DIR}/fibre_structure.hh
    ${SRC_DIR}/fibre_elongation_model.cc
    ${SRC_DIR}/fibre_elongation_model.hh
    ${SRC_DIR}/Fig_Object.cc
    ${SRC_DIR}/Fig_Object.hh
    ${SRC_DIR}/file.cc
    ${SRC_DIR}/file.hh
    ${SRC_DIR}/fixedwidth.hh
    ${SRC_DIR}/generic_synapse.cc
    ${SRC_DIR}/generic_synapse.hh
    ${SRC_DIR}/global.cc
    ${SRC_DIR}/global.hh
    ${SRC_DIR}/IFactivity.cc
    ${SRC_DIR}/IFactivity.hh
    ${SRC_DIR}/mtprng.cc
    ${SRC_DIR}/mtprng.hh
    ${SRC_DIR}/Network_Generated_Statistics.cc
    ${SRC_DIR}/Network_Generated_Statistics.hh
    ${SRC_DIR}/Network_Statistics.cc
    ${SRC_DIR}/Network_Statistics.hh

    ${SRC_DIR}/network.cc
    ${SRC_DIR}/network.hh
    ${SRC_DIR}/neurite_diameter_model.cc
    ${SRC_DIR}/neurite_diameter_model.hh
    ${SRC_DIR}/neuron.cc
    ${SRC_DIR}/neuron.hh
    # ${SRC_DIR}/nibr.cc
    # ${SRC_DIR}/nibr.hh
    ${SRC_DIR}/prepost_structure.cc
    ${SRC_DIR}/prepost_structure.hh
    ${SRC_DIR}/ptest.c

    ${SRC_DIR}/regex-gnu.c
    ${SRC_DIR}/regex-gnu.h
    ${SRC_DIR}/rx.h

    ${SRC_DIR}/Results.cc
    ${SRC_DIR}/Results.hh
    ${SRC_DIR}/Sampled_Output.cc
    ${SRC_DIR}/Sampled_Output.hh
    ${SRC_DIR}/slice.cc
    ${SRC_DIR}/slice.hh
    ${SRC_DIR}/Spatial_Presentation.cc
    ${SRC_DIR}/Spatial_Presentation.hh
    ${SRC_DIR}/Spatial_Segment_Subset.cc
    ${SRC_DIR}/Spatial_Segment_Subset.hh
    ${SRC_DIR}/spatial.cc
    ${SRC_DIR}/spatial.hh
    ${SRC_DIR}/state_storable.cc
    ${SRC_DIR}/state_storable.hh
    ${SRC_DIR}/StringList.cc
    ${SRC_DIR}/StringList.hh
    ${SRC_DIR}/synapse_formation_model.cc
    ${SRC_DIR}/synapse_formation_model.hh
    ${SRC_DIR}/synapse_structure.cc
    ${SRC_DIR}/synapse_structure.hh
    ${SRC_DIR}/synapse.cc
    ${SRC_DIR}/synapse.hh
    ${SRC_DIR}/turning_models.cc
    ${SRC_DIR}/turning_models.hh
    ${SRC_DIR}/Txt_Object.cc
    ${SRC_DIR}/Txt_Object.hh

    ${SRC_DIR}/vector.c
    ${SRC_DIR}/vector.h
    ${SRC_DIR}/Vector2D.cc
    ${SRC_DIR}/Vector2D.hh
    ${SRC_DIR}/Vector3D.cc
    ${SRC_DIR}/Vector3D.hh

    ${SRC_DIR}/VRML_Object.cc
    ${SRC_DIR}/VRML_Object.hh

    ${SRC_DIR}/vtest.c

    ${SRC_DIR}/include/templates.hh


)


# Add Library for dil2al
add_library("NetmorphDIL2AL" ${DIL2AL_SOURCES})
target_include_directories("NetmorphDIL2AL" PUBLIC ${SRC_DIR}/dil2al)
target_compile_definitions("NetmorphDIL2AL" PRIVATE
    -D__DIL2AL__
    -DSTDC_HEADERS
    -D_CPP_REGEX
    -D_BIGREGEX_SAFE_MATCHES 
)

# Add Library for geometry
add_library("NetmorphGeometry" ${GEOMETRY_SOURCES})
target_include_directories("NetmorphGeometry" PUBLIC ${SRC_DIR}/geometry)
target_compile_definitions("NetmorphGeometry" PRIVATE
    -D__DIL2AL__
    -DSTDC_HEADERS
    -D_CPP_REGEX
    -D_BIGREGEX_SAFE_MATCHES 
)





# Create Main Executable
if (NOT ${BUILD_NETMORPH_AS_LIB})
    add_executable(${PROJECT_NAME} "Include/Embeddable.h" "nibr.hh" "nibr.cc" ${MAIN_SOURCES})
    set(EXECUTABLE_TYPE "Binary")
    target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

    target_link_libraries(${PROJECT_NAME} PRIVATE 
        "NetmorphDIL2AL"
        "NetmorphGeometry"
    )
    target_include_directories(${PROJECT_NAME} PUBLIC ${SRC_DIR})
    target_include_directories(${PROJECT_NAME} PRIVATE ${SRC_DIR}/include)

    target_compile_definitions(${PROJECT_NAME} PRIVATE



        -DVECTOR3D
        
        -D__DIL2AL__ 
        -D_ALT_REGEX_H=\"regex-gnu.h\"
        -D_USE_ALT_REGEX
        -D_CPP_REGEX
        -DEVALUATE_POSSIBLE_CONNECTION_PROFILING
        -DSYNAPTOGENESIS_AND_LOSS_INVENTORY
        -DSAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
        -DADM_PRIORITIZE_SPEED_OVER_SPACE
        -DTESTING_SIMPLE_ATTRACTION
        -DPRECISE_RNG
        -DINCLUDE_SIMPLE_FIBER_DIAMETER
        -DUSE_RCFILE
        -DUSE_RCFILE
        -DENABLE_FIXED_STEP_SIMULATION
        -DFIXED_STEP_SIMULATION_START_AT_DT
        -DTESTING_SPIKING
        -DCHECK_XFIG_RANGES
        -D_BIGREGEX_SAFE_MATCHES
        -DFIG2DEV=\"fig2dev\"
        
    )
else()
    add_library(${PROJECT_NAME} "Include/Embeddable.h" "Include/Netmorph.h" "Include/Netmorph.cpp" ${MAIN_SOURCES})

    target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

    target_link_libraries(${PROJECT_NAME} PUBLIC 
        "NetmorphDIL2AL"
        "NetmorphGeometry"
    )
    target_include_directories(${PROJECT_NAME} PUBLIC ${SRC_DIR})
    target_include_directories(${PROJECT_NAME} PUBLIC ${SRC_DIR}/include)

    target_compile_definitions(${PROJECT_NAME} PUBLIC



    -DVECTOR3D
    
    -D__DIL2AL__ 
    -D_ALT_REGEX_H=\"regex-gnu.h\"
    -D_USE_ALT_REGEX
    -D_CPP_REGEX
    -DEVALUATE_POSSIBLE_CONNECTION_PROFILING
    -DSYNAPTOGENESIS_AND_LOSS_INVENTORY
    -DSAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
    -DADM_PRIORITIZE_SPEED_OVER_SPACE
    -DTESTING_SIMPLE_ATTRACTION
    -DPRECISE_RNG
    -DINCLUDE_SIMPLE_FIBER_DIAMETER
    -DUSE_RCFILE
    -DUSE_RCFILE
    -DENABLE_FIXED_STEP_SIMULATION
    -DFIXED_STEP_SIMULATION_START_AT_DT
    -DTESTING_SPIKING
    -DCHECK_XFIG_RANGES
    -D_BIGREGEX_SAFE_MATCHES
    -DFIG2DEV=\"fig2dev\"
    )
   
endif()