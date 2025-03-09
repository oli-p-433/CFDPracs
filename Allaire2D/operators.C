#include "operators.H"
#include <cmath>

// ------------------- Overloading operators ------------------ //
template<size_t N>
std::array<double, N> operator-(const std::array<double, N>& lhs, const std::array<double, N>& rhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}

template<size_t N>
std::array<double, N> operator+(const std::array<double, N>& lhs, const std::array<double, N>& rhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}

template<size_t N>
std::array<double, N> operator*(const double scalar,const std::array<double, N>& lhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = scalar*lhs[i];
    }
    return result;
}

template<size_t N>
std::array<double, N> operator*(const std::array<double, N>& rhs,const std::array<double, N>& lhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = rhs[i]*lhs[i];
    }
    return result;
}

template<size_t N>
std::array<double, N> elementDivide(const std::array<double, N>& lhs, const std::array<double, N>& rhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        if (rhs[i] == lhs[i]){
            result[i] = 1;
        }
        else if (rhs[i] == 0){
            result[i] = INFINITY;
        }
        else{
            result[i] = lhs[i] / rhs[i];
        }
    }
    return result;
}

template<size_t N>
std::array<double, N> operator/(const std::array<double, N>& lhs, const double scalar) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i]/scalar;
    }
    return result;
}