if(${CERES_FOUND})

    # Add inverse transfrom library
    add_library(inverse_compund_cf SHARED inverse_compound_cost_function.cpp)
    target_link_librarieS(inverse_compund_cf ceres glog)

    # Add inverse transfrom library
    add_library(compund_cf SHARED compound_cost_function.cpp)
    target_link_librarieS(compund_cf ceres glog)

    # Add quaternion from rpy library
    add_library(quaternion_from_rpy_cf SHARED quaternion_from_rpy_cost_function.cpp)
    target_link_librarieS(quaternion_from_rpy_cf ceres glog)

    # Add quaternion from rpy library
    add_library(quaternion_angle_cf SHARED quaternion_angle_cost_function.cpp)
    target_link_librarieS(quaternion_angle_cf ceres glog)

    # Add quaternion from rpy library
    add_library(quaternion_distance_cf SHARED quaternion_distance_cost_function.cpp)
    target_link_librarieS(quaternion_distance_cf ceres glog)

    # Add test
    add_executable(test_jacobians test_jacobians.cpp)
    target_link_libraries(test_jacobians compund_cf inverse_compund_cf quaternion_from_rpy_cf quaternion_angle_cf quaternion_distance_cf ceres glog )
endif()
