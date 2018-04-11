if(${CERES_FOUND})

    # Add inverse transfrom library
    add_library(inverse_compund_cf SHARED inverse_compound_cost_function.cpp)
    target_link_librarieS(inverse_compund_cf ceres glog)

    # Add inverse transfrom library
    add_library(compund_cf SHARED compound_cost_function.cpp)
    target_link_librarieS(compund_cf ceres glog)

    # Add test
    add_executable(test_jacobians test_jacobians.cpp)
    target_link_libraries(test_jacobians compund_cf inverse_compund_cf ceres glog )
endif()
