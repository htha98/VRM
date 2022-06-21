#include <random>
#include "kinematics.h"
//#define SHOW_ITERATION

/**
 * \brief Forward kinematics - DH matrices method
 * \details Target position is stored in the last column of the DH matrix of the target. \f$T=\prod_(i=1)^n L_i\f$.
 */
dh_matrix forward_kinematics(dh_matrix *p_dh_matrices, size_t num_of_matrices)
{
    dh_matrix target = p_dh_matrices[0];
    for(size_t i = 1; i < num_of_matrices; ++i)
        target *= p_dh_matrices[i];
    return target;
}

/**
 * \brief Inverse kinematics - FABRIK algorithm (https://doi.org/10.1016/j.gmod.2011.05.003)
 * \details FABRIK algorithm is heuristic method for fast iterative solving inverse kinematics problem. This function
 * (unlike the other one for calculating inverse kinematics problem) works with links of the same length. There is
 * also randomizer for faster iteration - it is off by default.
 */
point *inverse_kinematics(point start, point target, double link_len, size_t num_of_links,
                          bool set_random_coordinates=false, double tolerance=0.25, size_t max_iter=100)
{
#ifdef SHOW_ITERATION
    size_t iter_counter = 0;
#endif
    size_t num_of_joints = num_of_links + 1;
    // start with huge error
    double error = 1000 * tolerance;
    point direction = target - start;
    double distance = direction.magnitude();
    // check whether it is possible to reach the target
    if(distance > link_len * double(num_of_links))
        throw std::out_of_range("Target cannot be reached.");
    // work with joints (as points) instead of links
    auto *joints = new point[num_of_joints];
    // check whether randomizer is on, if so, generate random positions of joints (points)
    if(set_random_coordinates)
    {
        point min_corner((start.x < target.x ? start.x : target.x),
                         (start.y < target.y ? start.y : target.y),
                         (start.z < target.z ? start.z : target.z));
        point max_corner((start.x > target.x ? start.x : target.x),
                         (start.y > target.y ? start.y : target.y),
                         (start.z > target.z ? start.z : target.z));
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis_x(min_corner.x, max_corner.x);
        std::uniform_real_distribution<> dis_y(min_corner.y, max_corner.y);
        std::uniform_real_distribution<> dis_z(min_corner.z, max_corner.z);
        for (size_t i = 0; i < num_of_joints; ++i)
            joints[i] = {dis_x(gen), dis_y(gen), dis_z(gen)};
    }

    // temporary variables
    point tmp_vector;
    double tmp_magnitude;

    // calculate until either max. iteration reached or error is smaller than threshold value
    for(size_t iteration = 0; iteration < max_iter && error > tolerance; ++iteration)
    {
        // set target as the last point
        joints[num_of_joints - 1] = target;
        // backwards reaching inverse kinematics
        for (size_t i = 0; i < num_of_joints - 1; ++i)
        {
            // find the direction
            tmp_vector = joints[(num_of_joints - 1) - i - 1] - joints[(num_of_joints - 1) - i];
            // calculate the expand value
            tmp_magnitude = tmp_vector.magnitude();
            tmp_vector *= link_len / tmp_magnitude;
            // calculate the coordinates
            joints[(num_of_joints - 1) - i - 1] = joints[(num_of_joints - 1) - i] + tmp_vector;
        }

        // set base as the first point
        joints[0] = start;
        // forward reaching inverse kinematics
        for (size_t i = 0; i < num_of_joints - 1; ++i)
        {
            // find the direction
            tmp_vector = joints[i + 1] - joints[i];
            // calculate the expand value
            tmp_magnitude = tmp_vector.magnitude();
            tmp_vector *= link_len / tmp_magnitude;
            // calculate the coordinates
            joints[i + 1] = joints[i] + tmp_vector;
        }
        // calculate the error
        error = (joints[num_of_joints - 1] - target).magnitude();
#ifdef SHOW_ITERATION
        iter_counter = iteration;
#endif
    }
#ifdef SHOW_ITERATION
    std::cout << "\tNumber of iterations: " << iter_counter << std::endl;
#endif
    return joints;
}

/**
 * \brief Inverse kinematics - FABRIK algorithm (https://doi.org/10.1016/j.gmod.2011.05.003)
 * \details FABRIK algorithm is heuristic method for fast iterative solving inverse kinematics problem. This function
 * (unlike the other one for calculating inverse kinematics problem) works with links of different lengths. There is
 * also randomizer for faster iteration - it is off by default.
 */
point *inverse_kinematics(point start, point target, const double *link_len, size_t num_of_links,
                          bool set_random_coordinates=false, double tolerance=0.25, size_t max_iter=1000)
{
#ifdef SHOW_ITERATION
    size_t iter_counter = 0;
#endif
    size_t num_of_joints = num_of_links + 1;
    // start with huge error
    double error = 1000 * tolerance;
    point direction = target - start;
    double distance = direction.magnitude();
    double links_len = 0;
    // calculate the sum length of links
    for (size_t i = 0; i < num_of_links; ++i)
        links_len += link_len[i];
    // check whether it is possible to reach the target
    if(distance > links_len)
        throw std::out_of_range("Target cannot be reached.");
    // work with joints (as points) instead of links
    auto *joints = new point[num_of_joints];
    // check whether randomizer is on, if so, generate random positions of joints (points)
    if(set_random_coordinates)
    {
        point min_corner((start.x < target.x ? start.x : target.x),
                         (start.y < target.y ? start.y : target.y),
                         (start.z < target.z ? start.z : target.z));
        point max_corner((start.x > target.x ? start.x : target.x),
                         (start.y > target.y ? start.y : target.y),
                         (start.z > target.z ? start.z : target.z));
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis_x(min_corner.x, max_corner.x);
        std::uniform_real_distribution<> dis_y(min_corner.y, max_corner.y);
        std::uniform_real_distribution<> dis_z(min_corner.z, max_corner.z);
        for (size_t i = 0; i < num_of_joints; ++i)
            joints[i] = {dis_x(gen), dis_y(gen), dis_z(gen)};
    }
    // temporary variables
    point tmp_vector;
    double tmp_magnitude;

    // calculate until either max. iteration reached or error is smaller than threshold value
    for(size_t iteration = 0; iteration < max_iter && error > tolerance; ++iteration)
    {
        // set target as the last point
        joints[num_of_joints - 1] = target;
        // backwards reaching inverse kinematics
        for (size_t i = 0; i < num_of_joints - 1; ++i)
        {
            // find the direction
            tmp_vector = joints[(num_of_joints - 1) - i - 1] - joints[(num_of_joints - 1) - i];
            // calculate the expand value
            tmp_magnitude = tmp_vector.magnitude();
            tmp_vector *= link_len[num_of_links - 1 - i] / tmp_magnitude;
            // calculate the coordinates
            joints[(num_of_joints - 1) - i - 1] = joints[(num_of_joints - 1) - i] + tmp_vector;
        }

        // set base as the first point
        joints[0] = start;
        // forward reaching inverse kinematics
        for (size_t i = 0; i < num_of_joints - 1; ++i)
        {
            // find the direction
            tmp_vector = joints[i + 1] - joints[i];
            // calculate the expand value
            tmp_magnitude = tmp_vector.magnitude();
            tmp_vector *= link_len[i] / tmp_magnitude;
            // calculate the coordinates
            joints[i + 1] = joints[i] + tmp_vector;
        }
        // calculate the error
        error = (joints[num_of_joints - 1] - target).magnitude();
#ifdef SHOW_ITERATION
        iter_counter = iteration;
#endif
    }
#ifdef SHOW_ITERATION
    std::cout << "\tNumber of iterations: " << iter_counter << std::endl;
#endif
    return joints;
}

int main()
{
    // conversion ratio for converting degrees to radians
    constexpr double conv_ratio = M_PI / 180;
    // number of links
    const size_t n = 6;
    std::cout << "FORWARD KINEMATICS:" << std::endl;
    dh_matrix forward_matrices[n];
    // define parameters for each link
    dh_params forward_params[n] = {dh_params(290.0, 0.0, 0.0, -90.0 * conv_ratio),
                                   dh_params(0.0, -90.0 * conv_ratio, 270.0, 0.0),
                                   dh_params(0.0, 180.0 * conv_ratio, -70.0, 90.0 * conv_ratio),
                                   dh_params(302.0, 0.0, 0.0, -90.0 * conv_ratio),
                                   dh_params(0.0, 0.0, 0.0, 90.0 * conv_ratio),
                                   dh_params(72.0, 0.0, 0.0, 0.0)};
    // pass the parameters to the dh_matrix and calculate the matrices
    for (size_t i = 0; i < n; ++i)
        forward_matrices[i].set_params(forward_params[i]);
    // calculate DH matrix for the target
    dh_matrix target_fk = forward_kinematics(forward_matrices, n);
    double *target_matrix_fk = target_fk.get_matrix();
    // coordinates are x = T_{1,4}, y = T_{2,4} and z = T_{3,4}
    point target_coordinates_fk(round(target_matrix_fk[0 * 4 + 3]),
                                round(target_matrix_fk[1 * 4 + 3]),
                                round(target_matrix_fk[2 * 4 + 3]));
    std::cout << "Target DH matrix:" << std::endl;
    // print the matrix in a correct matrix format
    for (size_t i = 0; i < 4 * 4; ++i)
    {
        if(i == 4 || i == 2 * 4 || i == 3 * 4)
            std::cout << std::endl;
        std::cout << round(target_matrix_fk[i]) << "\t";
    }
    // print target coordinates
    std::cout << "\nTarget coordinates: " << target_coordinates_fk << std::endl;

    std::cout << "--------------------------------------" << std::endl;

    std::cout << "INVERSE KINEMATICS:" << std::endl;
    point target_coordinates_ik = target_coordinates_fk;
    // set_randomness for generating random points
    bool set_randomness = true;
//    bool set_randomness = false;
    std::cout << "- same length of links:" << std::endl;
    point base(0.0, 0.0, 0.0);
    point *joints = inverse_kinematics(base, target_coordinates_ik, 150, n, set_randomness);
    for (size_t i = 0; i < n + 1; ++i)
        if(i == 0)
            std::cout << "\tBase: " << joints[i] << std::endl;
        else if(i == n)
            std::cout << "\tTarget: " << joints[i] << std::endl;
        else
            std::cout << "\tJoint" << i + 1 << ": " << joints[i] << std::endl;

    std::cout << "\n- variable lengths of links:" << std::endl;
    double links_len[n] = {120, 150, 150, 120, 200, 100};
    joints = inverse_kinematics(base, target_coordinates_ik, links_len, n, set_randomness);
    for (size_t i = 0; i < n + 1; ++i)
        if(i == 0)
            std::cout << "\tBase: " << joints[i] << std::endl;
        else if(i == n)
            std::cout << "\tTarget: " << joints[i] << std::endl;
        else
            std::cout << "\tJoint" << i + 1 << ": " << joints[i] << std::endl;


    return 0;
}
