#include <iostream>

#include "../include/linear_algebra_system.hpp"

int main()
{
    namespace cramers_method = math::linear_algebra_system::cramers_method;

    math::linear_algebra_system::matrix_3x3<int64_t> matrix {};
    
    for (uint8_t i {0}; i < 3; i++) {
        std::cout << "Enter 3 values of " << i + 1 << " row separated by space: ";
        std::cin >> matrix[i][0] >> matrix[i][1] >> matrix[i][2];
    }

    std::cout << "Enter value of b1, b2 and b3 separated by space: ";
    int64_t b1 {}, b2 {}, b3 {};
    std::cin >> b1 >> b2 >> b3;

    const auto [X, Y, Z] {cramers_method::compute_all(matrix, {b1, b2, b3})};

    std::cout << "Answer: (X=" << X << "; Y=" << Y << "; Z=" << Z << ")\n";

    return 0;
}