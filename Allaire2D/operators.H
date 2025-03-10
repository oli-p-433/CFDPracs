#ifndef OPERATORS_H
#define OPERATORS_H

#include <array>
#include <iostream>
#include <cmath>

template<size_t N>
std::array<double, N> operator-(const std::array<double, N>& lhs, const std::array<double, N>& rhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] - rhs[i];
        if (std::isnan(result[i])) {
            std::cout << "LHS:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << lhs[j] << " ";
            }
            std::cout << std::endl;

            std::cout << "RHS:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << rhs[j] << " ";
            }
            std::cout << std::endl;

            std::cout << "Result:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << result[j] << " ";
            }
            std::cout << std::endl;

            throw std::runtime_error("element subtract error");
        }
    }
    return result;
}

template<size_t N>
std::array<double, N> operator+(const std::array<double, N>& lhs, const std::array<double, N>& rhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] + rhs[i];
        if (std::isnan(result[i])) {
            std::cout << "LHS:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << lhs[j] << " ";
            }
            std::cout << std::endl;

            std::cout << "RHS:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << rhs[j] << " ";
            }
            std::cout << std::endl;

            std::cout << "Result:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << result[j] << " ";
            }
            std::cout << std::endl;

            throw std::runtime_error("element add error");
        }
    }
    return result;
}

template<size_t N>
std::array<double, N> operator*(const double scalar, const std::array<double, N>& lhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = scalar * lhs[i];
        if (std::isnan(result[i])) {
            std::cout << "scalar " << scalar << std::endl;
            std::cout << "LHS:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << lhs[j] << " ";
            }
            std::cout << std::endl;

            std::cout << "Result:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << result[j] << " ";
            }
            std::cout << std::endl;

            throw std::runtime_error("scalar multiply error");
        }
    }
    return result;
}

template<size_t N>
std::array<double, N> operator*(const std::array<double, N>& rhs, const std::array<double, N>& lhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = rhs[i] * lhs[i];
        if (std::isnan(result[i])) {
            std::cout << "LHS:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << lhs[j] << " ";
            }
            std::cout << std::endl;

            std::cout << "RHS:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << rhs[j] << " ";
            }
            std::cout << std::endl;

            std::cout << "Result:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << result[j] << " ";
            }
            std::cout << std::endl;

            throw std::runtime_error("element multiply error");
        }
    }
    return result;
}

template<size_t N>
std::array<double, N> elementDivide(const std::array<double, N>& lhs, const std::array<double, N>& rhs) {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        if (rhs[i] == lhs[i]) {
            result[i] = 1;
        } else if (rhs[i] == 0) {
            result[i] = INFINITY;
        } else {
            result[i] = lhs[i] / rhs[i];
        }
        if (std::isnan(result[i])) {
            std::cout << "LHS:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << lhs[j] << " ";
            }
            std::cout << std::endl;

            std::cout << "RHS:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << rhs[j] << " ";
            }
            std::cout << std::endl;

            std::cout << "Result:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << result[j] << " ";
            }
            std::cout << std::endl;

            throw std::runtime_error("element divide error");
        }
    }
    return result;
}

template<size_t N>
std::array<double, N> operator/(const std::array<double, N>& lhs, const double scalar) {
    if (scalar == 0) {
        throw std::runtime_error("Division by zero error");
    }

    std::array<double, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] / scalar;
        if (std::isnan(result[i])) {
            std::cout << "LHS:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << lhs[j] << " ";
            }
            std::cout << std::endl;

            std::cout << "Result:";
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << result[j] << " ";
            }
            std::cout << std::endl;

            throw std::runtime_error("Scalar divide error");
        }
    }
    return result;
}

template<typename T>
void resizeVector(std::vector<std::vector<T>>& matrix, size_t rows, size_t cols) {
    matrix.resize(rows);
    for(auto& row : matrix) {
        row.resize(cols);
    }
}

#endif // OPERATORS_H
