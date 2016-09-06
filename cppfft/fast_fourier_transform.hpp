#ifndef CPPFFT_FAST_FOURIER_TRANSFORM_HPP
#define CPPFFT_FAST_FOURIER_TRANSFORM_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#include "./detail/is_inputtable_iterator.hpp"
#include "./detail/is_outputtable_iterator.hpp"
#include "./detail/is_random_accessible_iterator.hpp"

namespace cppfft {

namespace detail {

template <typename InputIterator1, typename InputIterator2, typename ForwardIterator>
constexpr auto next_positional(InputIterator1 first, InputIterator2 last, ForwardIterator result)
{
    using value_type = typename std::iterator_traits<ForwardIterator>::value_type;
    for (;
        first != last && (*result += value_type(1)) >= *first;
        void(++first), *result++ = value_type(0));
    return result;
}

template <typename InputIterator, typename Size, typename ForwardIterator>
inline constexpr auto next_positional_n(InputIterator first, Size size, ForwardIterator result)
{
    return cppfft::detail::next_positional(first, first + size, result);
}

template <typename DifferenceType, typename RandomAccessIterator1, typename RandomAccessIterator2>
auto replace(
    std::vector<DifferenceType> const& radices,
    RandomAccessIterator1 first,
    DifferenceType size,
    RandomAccessIterator2 result)
    -> std::enable_if_t<
        cppfft::detail::is_random_accessible_iterator_v<RandomAccessIterator1>,
        RandomAccessIterator2>
{
    using difference_type = DifferenceType;
    auto const& length = radices.back();
    auto const stride = size / length;

    auto coefficient = std::vector<difference_type>(radices.size() - 1u);
    coefficient.back() = difference_type{1};

    for (auto i = std::size_t{1}; i < coefficient.size(); ++i)
    {
        auto const index = coefficient.size() - 1u - i;
        coefficient.at(index) = coefficient.at(index + 1u) * radices.at(index + 1u);
    }

    auto positions = std::vector<difference_type>(radices.size() - 1u, difference_type{0});

    for (auto i = difference_type{0}; i < stride; void(++i),
        cppfft::detail::next_positional_n(
            radices.cbegin(), positions.size(), positions.begin()))
    {
        auto const index = length * std::inner_product(
            positions.cbegin(), positions.cend(), coefficient.cbegin(), difference_type{0});

        for (auto j = difference_type{0}; j < length; ++j)
        {
            result[index + j] = first[i + stride * j];
        }
    }

    return result + size;
}

template <typename DifferenceType, typename InputIterator, typename RandomAccessIterator>
auto replace(
    std::vector<DifferenceType> const& radices,
    InputIterator first,
    DifferenceType size,
    RandomAccessIterator result)
    -> std::enable_if_t<
        !cppfft::detail::is_random_accessible_iterator_v<InputIterator>
            && cppfft::detail::is_inputtable_iterator_v<InputIterator>,
        RandomAccessIterator>
{
    using difference_type = DifferenceType;
    auto const length = radices.back();
    auto const stride = size / length;

    std::vector<difference_type> coefficient(radices.size() - 1u);
    coefficient.back() = difference_type{1};

    for (auto i = std::size_t{1}; i < coefficient.size(); ++i)
    {
        auto const index = coefficient.size() - 1u - i;
        coefficient.at(index) = coefficient.at(index + 1u) * radices.at(index + 1u);
    }

    auto positions = std::vector<difference_type>(radices.size() - 1u, difference_type{0});

    for (auto i = difference_type{0}; i < length; ++i)
    {
        for (auto j = difference_type{0}; j < stride;
            void(++j), void(++first),
            cppfft::detail::next_positional_n(
                radices.cbegin(), positions.size(), positions.begin()))
        {
            auto const index = length * std::inner_product(
                positions.cbegin(), positions.cend(), coefficient.cbegin(), difference_type{0});
            result[index + i] = *first;
        }
    }

    return result + size;
}

} // namespace detail

constexpr struct forward_t {} forward{};
constexpr struct inverse_t {} inverse{};

template <typename ComplexType, typename DifferenceType = std::ptrdiff_t>
class fast_fourier_transform
{
public:
    using value_type = ComplexType;
    using element_type = typename value_type::value_type;
    using difference_type = DifferenceType;

private:
    std::vector<value_type> twiddles;
    std::vector<difference_type> radices;
    difference_type sequence_size;

public:
    fast_fourier_transform() = delete;
    fast_fourier_transform(fast_fourier_transform const&) = default;
    fast_fourier_transform(fast_fourier_transform&&) = default;

    explicit fast_fourier_transform(difference_type n)
        : sequence_size(n)
    {
        if (n < 0)
        {
            throw std::length_error(
                "fast_fourier_transform::fast_fourier_transform: size must be non-negative\n");
        }

        twiddles.resize(sequence_size);

        using std::acos;
        auto const k = element_type(-2. * acos(element_type(-1.)) / element_type(n));

        for (auto i = std::size_t{0u}; i < sequence_size; ++i)
        {
            using std::exp;
            twiddles.at(i) = exp(value_type(element_type(0.), element_type(i) * k));
        }

        for (auto i = difference_type{4}; i * i <= n; )
        {
            if (n % i != 0)
            {
                i = i == 2 ? 3 : i == 4 ? 2 : i + 2;
                continue;
            }

            radices.push_back(i);
            n /= i;
        }

        radices.push_back(n);
    }

    ~fast_fourier_transform() = default;

    auto operator=(fast_fourier_transform const&) & -> fast_fourier_transform& = default;
    auto operator=(fast_fourier_transform&&) & -> fast_fourier_transform& = default;

    auto size() const noexcept -> difference_type const&
    {
        return sequence_size;
    }

    template <typename TransformType, typename InputIterator, typename OutputIterator>
    auto operator()(TransformType&& type, InputIterator first, OutputIterator result) const
        -> std::enable_if_t<
            (std::is_same<std::decay_t<TransformType>, cppfft::forward_t>::value
                || std::is_same<std::decay_t<TransformType>, cppfft::inverse_t>::value)
                && cppfft::detail::is_inputtable_iterator_v<InputIterator>
                && cppfft::detail::is_outputtable_iterator_v<OutputIterator>,
            OutputIterator>
    {
        auto const is_inverse
            = std::is_same<std::decay_t<TransformType>, cppfft::inverse_t>::value;

        auto buffer = std::vector<value_type>(sequence_size);

        cppfft::detail::replace(radices, first, sequence_size, buffer.begin());

        auto remainders = std::vector<difference_type>(radices.size());
        auto strides = std::vector<difference_type>(radices.size());

        remainders.front() = sequence_size / radices.front();
        strides.front() = difference_type{1};

        for (auto i = std::size_t{0u}; i + 1 < radices.size(); ++i)
        {
            remainders.at(i + 1u) = remainders.at(i) / radices.at(i + 1u);
            strides.at(i + 1u) = strides.at(i) * radices.at(i);
        }

        for (auto i = std::size_t{0u}; i < radices.size(); ++i)
        {
            auto const index = radices.size() - i - 1u;
            auto const& stride = strides.at(index);
            auto const& radix = radices.at(index);
            auto const& remainder = remainders.at(index);

            for (auto iter = buffer.begin(); iter != buffer.end(); iter += radix * remainder)
            {
                butterfly(radix, is_inverse, remainder, stride, iter);
            }
        }

        return is_inverse
            ? std::transform(buffer.begin(), buffer.end(), result,
                [&](auto const& v) { return v / element_type(sequence_size); })
            : std::move(buffer.begin(), buffer.end(), result);
    }

    template <typename InputIterator, typename OutputIterator>
    auto operator()(InputIterator first, OutputIterator result) const
        -> std::enable_if_t<
            cppfft::detail::is_inputtable_iterator_v<InputIterator>
                && cppfft::detail::is_outputtable_iterator_v<OutputIterator>,
            OutputIterator>
    {
        return (*this)(cppfft::forward, first, result);
    }

    template <typename InputIterator, typename OutputIterator>
    auto forward(InputIterator first, OutputIterator result) const
        -> std::enable_if_t<
            cppfft::detail::is_inputtable_iterator_v<InputIterator>
                && cppfft::detail::is_outputtable_iterator_v<OutputIterator>,
            OutputIterator>
    {
        return (*this)(cppfft::forward, first, result);
    }

    template <typename InputIterator, typename OutputIterator>
    auto inverse(InputIterator first, OutputIterator result) const
        -> std::enable_if_t<
            cppfft::detail::is_inputtable_iterator_v<InputIterator>
                && cppfft::detail::is_outputtable_iterator_v<OutputIterator>,
            OutputIterator>
    {
        return (*this)(cppfft::inverse, first, result);
    }

private:
    template <typename T>
    static inline constexpr auto conjugate(bool const flag, T&& x) -> decltype(auto)
    {
        using std::conj;
        return flag ? conj(std::forward<T>(x)) : std::forward<T>(x);
    }

    auto butterfly(
        difference_type const& radix,
        bool const is_inverse,
        difference_type const& remainder,
        difference_type const& stride,
        typename std::vector<value_type>::iterator first) const
    {
        switch (radix)
        {
        case difference_type{2}: return butterfly2(is_inverse, remainder, stride, first);
        case difference_type{3}: return butterfly3(is_inverse, remainder, stride, first);
        case difference_type{4}: return butterfly4(is_inverse, remainder, stride, first);
        case difference_type{5}: return butterfly5(is_inverse, remainder, stride, first);
        default: break;
        }

        auto tmp = std::vector<value_type>(radix);

        for (auto i = difference_type{0}; i < remainder; ++i)
        {
            for (auto j = difference_type{0}; j < radix; ++j)
            {
                tmp.at(j) = first[i + j * remainder];
            }

            for (auto j = i; j < radix; j += remainder)
            {
                first[j] = value_type(0.);

                for (auto k = difference_type{0}; k < radix; ++k)
                {
                    first[j] += tmp.at(k)
                        * conjugate(is_inverse, twiddles.at((k * j * stride) % sequence_size));
                }
            }
        }
    }

    auto butterfly2(
        bool const is_inverse,
        difference_type const& remainder,
        difference_type const& stride,
        typename std::vector<value_type>::iterator first) const
    {
        for (auto i = difference_type{0}; i < remainder; ++i)
        {
            auto const t = first[remainder + i] * conjugate(is_inverse, twiddles.at(i * stride));
            first[remainder + i] = first[i] - t;
            first[i] += t;
        }
    }

    auto butterfly3(
        bool const is_inverse,
        difference_type const& remainder,
        difference_type const& stride,
        typename std::vector<value_type>::iterator first) const
    {
        using std::real;
        using std::imag;

        auto const t0 = imag(conjugate(is_inverse, twiddles.at(remainder * stride)));

        for (auto i = difference_type{0}; i < remainder; ++i)
        {
            auto const t1 = first[remainder + i] * conjugate(is_inverse, twiddles.at(i * stride));
            auto const t2
                = first[2 * remainder + i] * conjugate(is_inverse, twiddles.at(2 * i * stride));
            auto const t3 = t1 + t2;
            auto const t4 = (t1 - t2) * t0;
            first[remainder + i] = value_type(
                real(first[i]) - 0.5 * real(t3), imag(first[i]) - 0.5 * imag(t3));
            first[i] += t3;
            first[2 * remainder + i] = value_type(
                real(first[remainder + i]) + imag(t4), imag(first[remainder + i]) - real(t4));
            first[remainder + i] += value_type(-imag(t4), real(t4));
        }
    }

    auto butterfly4(
        bool const is_inverse,
        difference_type const& remainder,
        difference_type const& stride,
        typename std::vector<value_type>::iterator first) const
    {
        using std::real;
        using std::imag;

        for (auto i = difference_type{0}; i < remainder; ++i)
        {
            auto const t0
                = first[remainder + i] * conjugate(is_inverse, twiddles.at(i * stride));
            auto const t1
                = first[2 * remainder + i] * conjugate(is_inverse, twiddles.at(2 * i * stride));
            auto const t2
                = first[3 * remainder + i] * conjugate(is_inverse, twiddles.at(3 * i * stride));
            auto const t3 = first[i] - t1;
            first[i] += t1;
            auto const t4 = t0 + t2;
            auto t5 = t0 - t2;
            t5 = (is_inverse ? -1. : 1.) * value_type(imag(t5), -real(t5));
            first[2 * remainder + i] = first[i] - t4;
            first[i] += t4;
            first[remainder + i] = t3 + t5;
            first[3 * remainder + i] = t3 - t5;
        }
    }

    auto butterfly5(
        bool const is_inverse,
        difference_type const& remainder,
        difference_type const& stride,
        typename std::vector<value_type>::iterator first) const
    {
        using std::real;
        using std::imag;

        auto const t0 = conjugate(is_inverse, twiddles.at(remainder * stride));
        auto const t1 = conjugate(is_inverse, twiddles.at(2u * remainder * stride));

        for (auto i = difference_type{0}; i < remainder; ++i)
        {
            auto const t2 = first[i];
            auto const t3
                = first[remainder + i] * conjugate(is_inverse, twiddles.at(i * stride));
            auto const t4
                = first[2 * remainder + i] * conjugate(is_inverse, twiddles.at(2 * i * stride));
            auto const t5
                = first[3 * remainder + i] * conjugate(is_inverse, twiddles.at(3 * i * stride));
            auto const t6
                = first[4 * remainder + i] * conjugate(is_inverse, twiddles.at(4 * i * stride));
            auto const t7 = t3 + t6;
            auto const t8 = t3 - t6;
            auto const t9 = t4 + t5;
            auto const t10 = t4 - t5;
            first[i] = first[i] + t7 + t9;
            auto const t11 = t2 + value_type(
                real(t7) * real(t0) + real(t9) * real(t1),
                imag(t7) * real(t0) + imag(t9) * real(t1));
            auto const t12 = value_type(
                imag(t8) * imag(t0) + imag(t10) * imag(t1),
                -real(t8) * imag(t0) - real(t10) * imag(t1));
            first[remainder + i] = t11 - t12;
            first[4 * remainder + i] = t11 + t12;
            auto const t13 = t2 + value_type(
                real(t7) * real(t1) + real(t9) * real(t0),
                imag(t7) * real(t1) + imag(t9) * real(t0));
            auto const t14 = value_type(
                -imag(t8) * imag(t1) + imag(t10) * imag(t0),
                real(t8) * imag(t1) - real(t10) * imag(t0));
            first[2 * remainder + i] = t13 + t14;
            first[3 * remainder + i] = t13 - t14;
        }
    }
};

} // namespace cppfft

#endif // #ifndef CPPFFT_FAST_FOURIER_TRANSFORM_HPP
