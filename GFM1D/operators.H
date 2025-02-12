#include <array>

std::array<double, 3> operator-(const std::array<double, 3>&, const std::array<double, 3>&);

std::array<double, 3> operator+(const std::array<double, 3>&, const std::array<double, 3>&);

std::array<double, 3> operator*(const double,const std::array<double, 3>&);

std::array<double, 3> operator*(const std::array<double, 3>&,const std::array<double, 3>&);

std::array<double, 3> operator/(const std::array<double, 3>&, const double);

std::array<double, 3> elementDivide(const std::array<double, 3>&, const std::array<double, 3>&);