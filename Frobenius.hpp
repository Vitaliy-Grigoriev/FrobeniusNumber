#pragma once
#ifndef FROBENIUS_NUMBER_HPP
#define FROBENIUS_NUMBER_HPP

#include <vector>
#include <climits>
#include <iostream>
#include <algorithm>
#include <type_traits>


/**
 * @fn template <typename Type>
 * typename std::enable_if_t<std::is_arithmetic<Type>::value, Type>
 * constexpr GCD (const Type, const Type) noexcept;
 *
 * @brief Function that calculates the Greatest Common Divisor of two operands.
 * @tparam [in] Type - Template arithmetic value type for internal calculations.
 * @tparam [in] first - First operand.
 * @tparam [in] second - Second operand.
 * @return The Greatest Common Divisor result.
 */
template <typename Type>
typename std::enable_if_t<std::is_arithmetic<Type>::value, Type>
constexpr GCD (const Type first, const Type second) noexcept
{
    if (second == 0) { return first; }
    return GCD(second, first % second);
}

/**
 * @fn template <typename TypeIt>
 * typename std::iterator_traits<TypeIt>::value_type
 * constexpr CommonGCD (const TypeIt, const TypeIt) noexcept;
 *
 * @brief Function that calculates the Greatest Common Divisor of several operands.
 * @tparam [in] TypeIt - Type of iterator.
 * @tparam [in] first - An iterator to the first element of the container.
 * @tparam [in] last - An iterator to the end element of the container.
 * @return The Greatest Common Divisor result.
 */
template <typename TypeIt>
typename std::iterator_traits<TypeIt>::value_type
constexpr CommonGCD (const TypeIt first, const TypeIt last) noexcept
{
    const auto size = std::distance(first, last);
    if (size == 1) { return GCD(*first, *last); }
    const auto h = size / 2;
    return GCD(CommonGCD(first, first + h), CommonGCD(first + h, last));
}


/**
 * @fntemplate <typename Type>
 * typename std::enable_if_t<std::is_integral<Type>::value && std::is_signed<Type>::value, Type>
 * constexpr FrobeniusNumber (std::vector<Type> &) noexcept;
 *
 * @brief Function that calculates the Frobenius Number of the set of numbers.
 * @tparam [in] Type - Template arithmetic signed value type of input numbers.
 * @tparam [in] InitNumbersList - Input values for calculation.
 * @return The Frobenius Number of the set of numbers (Type::infinity - if an error occurred).
 */
template <typename Type>
typename std::enable_if_t<std::is_integral<Type>::value && std::is_signed<Type>::value, Type>
constexpr FrobeniusNumber (std::vector<Type>& InitNumbersList) noexcept
{
    if (CommonGCD(InitNumbersList.begin(), InitNumbersList.end() - 1) == Type(1))
    {
        std::sort(InitNumbersList.begin(), InitNumbersList.end());
        if (InitNumbersList[0] <= 0) { return std::numeric_limits<Type>::infinity(); }
        if (InitNumbersList[0] == 1) { return -1; }
        if (InitNumbersList.size() == 2) {
            return (InitNumbersList[0] * InitNumbersList[1] - InitNumbersList[0] - InitNumbersList[1]);
        }

        std::vector<Type> array(InitNumbersList[0], -1);
        array[0] = 0;

        for (std::size_t idx = 1; idx < InitNumbersList.size(); ++idx)
        {
            const Type d = GCD<Type>(InitNumbersList[0], InitNumbersList[idx]);
            for (Type r = 0; r < d; ++r)
            {
                Type n = -1;
                if (r == 0) { n = 0; }
                else
                {
                    Type q = r;
                    while (q < InitNumbersList[0])
                    {
                        if (array[q] != -1 && (array[q] < n || n == -1)) {
                            n = array[q];
                        }
                        q += d;
                    }
                }

                if (n != -1)
                {
                    for (std::size_t jdx = 0; jdx < InitNumbersList[0] / d; ++jdx)
                    {
                        n += InitNumbersList[idx];
                        const Type p = n % InitNumbersList[0];
                        if (array[p] != -1 && (array[p] < n || n == -1)) {
                            n = array[p];
                        }
                        array[p] = n;
                    }
                }
            }
        }

        Type max = 0;
        for (std::size_t idx = 0; idx < InitNumbersList[0]; ++idx) {
            if (array[idx] == -1 || array[idx] > max) {
                max = array[idx];
            }
        }
        return (max == -1) ? max : max - InitNumbersList[0];
    }
    return std::numeric_limits<Type>::infinity();
}

#endif  // FROBENIUS_NUMBER_HPP
