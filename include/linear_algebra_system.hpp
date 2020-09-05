#ifndef __LINEAR_ALGEBRA_SYSTEM_SOLE_CRAMERS_METHOD_HPP__
#define __LINEAR_ALGEBRA_SYSTEM_SOLE_CRAMERS_METHOD_HPP__

#include <cstdint>
#include <cstddef>

#include <array>
#include <tuple>

namespace math::linear_algebra_system
{
    template<typename MatrixType, size_t Rows, size_t Columns,
        std::enable_if_t<
            std::is_integral<MatrixType>::value,
            int64_t
        > = 0
    > using matrix = std::array<std::array<MatrixType, Rows>, Columns>;

    template<typename MatrixType,
        std::enable_if_t<
            std::is_integral<MatrixType>::value,
            int64_t
        > = 0
    > using matrix_3x3 = matrix<MatrixType, 3, 3>;

    namespace cramers_method
    {
        enum class DeltaType : uint8_t
        {
            X, Y, Z
        };

        template<typename MatrixType,
            std::enable_if_t<
                std::is_integral<MatrixType>::value,
                int64_t
            > = 0
        > constexpr auto
        compute_delta(const matrix_3x3<MatrixType> matrix)
        {
            if (!(matrix.size() >= 3 && matrix[0].size() >= 3)) {
                throw std::runtime_error("Matrix must be 3x3!");
            }

            return (
                ((matrix[0][0] * matrix[1][1] * matrix[2][2]) +
                (matrix[0][2] * matrix[2][1] * matrix[1][0]) +
                (matrix[0][1] * matrix[1][2] * matrix[2][0]))
                -
                ((matrix[0][2] * matrix[1][1] * matrix[2][0]) +
                (matrix[0][1] * matrix[1][0] * matrix[2][2]) +
                (matrix[0][0] * matrix[1][2] * matrix[2][1]))
            );
        }

        template<typename MatrixType,
            std::enable_if_t<
                std::is_integral<MatrixType>::value,
                int64_t
            > = 0
        > constexpr auto
        compute_delta(matrix_3x3<MatrixType> matrix, const DeltaType delta_type,
                        const std::tuple<MatrixType, MatrixType, MatrixType>& b_values)
        {
            if (!(matrix.size() >= 3 && matrix[0].size() >= 3)) {
                throw std::runtime_error("Matrix must be 3x3!");
            }

            switch (delta_type) {
            case DeltaType::X: std::tie(matrix[0][0], matrix[1][0], matrix[2][0]) = b_values; break;
            case DeltaType::Y: std::tie(matrix[0][1], matrix[1][1], matrix[2][1]) = b_values; break;
            case DeltaType::Z: std::tie(matrix[0][2], matrix[1][2], matrix[2][2]) = b_values; break;
            }

            return compute_delta(matrix);
        }

        template<typename MatrixType,
            std::enable_if_t<
                std::is_integral<MatrixType>::value,
                int64_t
            > = 0
        > inline constexpr auto
        compute_all(const matrix_3x3<MatrixType> matrix, const std::tuple<MatrixType, MatrixType, MatrixType>& b_values)
        {
            const auto delta {compute_delta(matrix)};
            return std::make_tuple(compute_delta(matrix, DeltaType::X, b_values) / delta,
                                    compute_delta(matrix, DeltaType::Y, b_values) / delta,
                                    compute_delta(matrix, DeltaType::Z, b_values) / delta);
        }
    } // namespace cramers_method
} // namespace math::linear_algebra_system

#endif // __LINEAR_ALGEBRA_SYSTEM_SOLE_CRAMERS_METHOD_HPP__