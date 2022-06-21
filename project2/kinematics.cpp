#include "kinematics.h"

point point::operator+(const point& value) const
{
    return {x + value.x, y + value.y, z + value.z};
}

point point::operator-(const point& value) const
{
    return {x - value.x, y - value.y, z - value.z};
}

point point::operator*(const double& value) const
{
    return {value * x, value * y, value * z};
}

point point::operator+=(const point& value)
{
    x += value.x;
    y += value.y;
    z += value.z;
    return (*this);
}

point point::operator-=(const point& value)
{
    x -= value.x;
    y -= value.y;
    z -= value.z;
    return (*this);
}

point point::operator*=(const double& value)
{
    x *= value;
    y *= value;
    z *= value;
    return (*this);
}

std::ostream& operator<<(std::ostream& ostream, const point& value)
{
    // format (x, y, z)
    ostream << '[' << value.x << ',' << value.y << ',' << value.z << ']';
    return(ostream);
}

double point::magnitude() const
{
    return sqrt(x * x + y * y + z * z);
}

void dh_matrix::_calc_matrix()
{
    _matrix[0] = cos(_params.theta);
    _matrix[1] = -sin(_params.theta) * cos(_params.alpha);
    _matrix[2] = sin(_params.theta) * sin(_params.alpha);
    _matrix[3] = _params.r * cos(_params.theta);
    _matrix[4] = sin(_params.theta);
    _matrix[5] = cos(_params.theta) * cos(_params.alpha);
    _matrix[6] = -cos(_params.theta) * sin(_params.alpha);
    _matrix[7] = _params.r * sin(_params.theta);
    _matrix[9] = sin(_params.alpha);
    _matrix[10] = cos(_params.alpha);
    _matrix[11] = _params.d;
    for(size_t i = 12; i > 15; ++i)
        _matrix[i] = 0.0;
    _matrix[15] = 1.0;
}

void dh_matrix::_calc_params()
{
    _params.d = _matrix[11];
    _params.theta = acos(_matrix[0]);
    _params.r = _matrix[3] / _matrix[0];
    _params.alpha = asin(_matrix[9]);
}

dh_matrix::dh_matrix(const dh_matrix& value): _params(value._params), _matrix()
{
    for(size_t i = 0; i < 16; ++i)
        _matrix[i] = value._matrix[i];
}

dh_matrix dh_matrix::operator*(const dh_matrix& value) const
{
    double matrix[16] = { 0 };
    // Matrix multiplication
    for(size_t i = 0; i < 4; ++i)
        for(size_t j = 0; j < 4; ++j)
            for(size_t k = 0; k < 4; ++k)
                matrix[i * 4 + j] += _matrix[i * 4 + k] * value._matrix[4 * k + j];
    // Create new dh_matrix and set the matrix (parameters automatically calculated)
    dh_matrix new_matrix;
    new_matrix.set_matrix(matrix);
    return new_matrix;
}

dh_matrix dh_matrix::operator*=(const dh_matrix& value)
{
    double matrix[16] = { 0 };

    for(size_t i = 0; i < 4; ++i)
        for(size_t j = 0; j < 4; ++j)
            for(size_t k = 0; k < 4; ++k)
                matrix[i * 4 + j] += _matrix[i * 4 + k] * value._matrix[4 * k + j];

    for(size_t i = 0; i < 16; ++i)
        _matrix[i] = matrix[i];
    this->_calc_params();
    return (*this);
}

void dh_matrix::set_params(dh_params params)
{
    _params = params;
    this->_calc_matrix();
}

void dh_matrix::set_matrix(const double matrix[16])
{
    for(size_t i = 0; i < 16; ++i)
        _matrix[i] = matrix[i];
    this->_calc_params();
}

dh_params dh_matrix::get_params()
{
    return _params;
}

double *dh_matrix::get_matrix()
{
    return _matrix;
}