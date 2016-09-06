#ifndef CPPFFT_IFFT_N_AS_HPP
#define CPPFFT_IFFT_N_AS_HPP

#include <iterator>
#include <type_traits>
#include "./detail/is_inputtable_iterator.hpp"
#include "./fast_fourier_transform.hpp"

namespace cppfft {

template <
    typename ComplexType,
    typename InputIterator,
    typename DifferenceType,
    typename OutputIterator>
inline auto ifft_n_as(InputIterator first, DifferenceType size, OutputIterator result)
    -> std::enable_if_t<
        cppfft::detail::is_inputtable_iterator_v<InputIterator>,
        OutputIterator>
{
    return cppfft::fast_fourier_transform<ComplexType, DifferenceType>(size)
        .inverse(first, result);
}

} // namespace cppfft

#endif // #ifndef CPPFFT_IFFT_N_AS_HPP
