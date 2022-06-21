#ifndef KINEMATICS_H
#define KINEMATICS_H
#include <cmath>
#include <iostream>

using std::sin;
using std::asin;
using std::cos;
using std::acos;
using std::sqrt;

/**
 * \brief point structure
 * \details Structure `point` defined by 3 values of the coordinate `x`, `y` and `z`. Some methods (operations) are also defined.
 */
struct point
{
    double x, y, z;
    /**
     * \brief Implicit constructor
     *  \details Initialize all values to 0.
     */
    point(): x(0), y(0), z(0) {}

    /**
     * \brief Conversion constructor
     * \details Copies double values to elements of the point.
     */
    point(double a_x, double a_y, double a_z): x(a_x), y(a_y), z(a_z) {}

    /**
     * \brief Operator +
     * \details Addition of two points is defined as addition of their elements.
     */
    point operator+(const point& value) const;

    /**
     * \brief Operator -
     * \details Subtraction of two points is defined as subtraction of their elements.
     */
    point operator-(const point& value) const;

    /**
     * \brief Operator * - multiplication of point by a single number
     * \details Multiplication of point `p` by a number `a` is defined as multiplying every element of point `p` by number `a`.
     */
    point operator*(const double& value) const;

    /**
     * \brief Operator +=
     * \details Almost same as operator `+` except the result of addition is saved to the variable `a` for `a += b`.
     */
    point operator+=(const point& value);

    /**
     * \brief Operator -=
     * \details Almost same as operator `-` except the result of subtraction is saved to the variable `a` for `a -= b`.
     */
    point operator-=(const point& value);

    /**
     * \brief Operator *= - multiplication of point by a single number
     * \details Defined as an enlargement/expansion of the point.
     */
    point operator*=(const double& value);

    /**
     * \brief Operator output to the stream <<
     * \details Defines the output to the stream as `(x,y,z)`.
     */
    friend std::ostream& operator<<(std::ostream& ostream, const point& value);

    /**
     * \brief Method for calculating magnitude of a vector
     * \details Magnitude \f$|v|\f$ defined as Euclidean distance: \f$|v|=\sqrt{x^2+y^2+z^2}\f$.
     */
    double magnitude() const;
};

/**
 * \brief Structure for storing Denavit-Hartenberg parameters
 * \details Structure `dh_params` contains parameters `d`, `theta`, `r` and `alpha`. Two constructors are defined.
 */
struct dh_params
{
    double d, theta, r, alpha;

    /**
     * \brief Implicit constructor
     * \details Initialize all values to 0.
     */
    dh_params(): d(0), theta(0), r(0), alpha(0) {}

    /**
     * \brief Conversion constructor
     * \details Copies all values to elements of the structure dh_params.
     */
    dh_params(double a_d, double a_theta, double a_r, double a_alpha): d(a_d), theta(a_theta), r(a_r), alpha(a_alpha) {}
};

/**
 * \brief Denevit-Hartenberg matrix
 * \details Matrix defined as array of n * m elements. Since n = m = 4, size of array is 16. For more details
 * about DH matrix: https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters#Denavit%E2%80%93Hartenberg_matrix
 */
class dh_matrix
{
    dh_params _params;
    double _matrix[16];

    /**
     * \brief Matrix calculator
     * \details Calculates matrix from parameters.
     */
    void _calc_matrix();

    /**
     * \brief Parameters calculator
     * \details Calculates parameters from matrix.
     */
    void _calc_params();
public:

    /**
     * \brief Implicit constructor
     * \details Initialize all values to 0.
     */
    dh_matrix(): _params(), _matrix() {}

    /**
     * \brief Copy constructor
     * \details Copies all values to the new dh_matrix.
     */
    dh_matrix(const dh_matrix& value);

    /**
     * \brief Operator *
     * \details A new dh_matrix contains matrix, which is result of matrix multiplication and then parameters are calculated.
     */
    dh_matrix operator*(const dh_matrix& value) const;

    /**
     * \brief Operator *=
     * \details Similar to operator *, however new dh_matrix is stored in place of the one used for the operation.
     */
    dh_matrix operator*=(const dh_matrix& value);

    /**
     * \brief Parameters setter
     * \details Sets parameters of dh_matrix based on input argument.
     */
    void set_params(dh_params params);

    /**
     * \brief Matrix setter
     * \details Sets matrix of dh_matrix based on input argument.
     */
    void set_matrix(const double matrix[16]);

    /**
     * \brief Parameters getter
     * \details Returns parameters of dh_matrix.
     */
    dh_params get_params();

    /**
     * \brief Matrix getter
     * \details Returns matrix of dh_matrix.
     */
    double* get_matrix();
};
#endif //KINEMATICS_H
