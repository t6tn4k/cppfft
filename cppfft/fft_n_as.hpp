#ifndef CPPFFT_FFT_N_AS_HPP
#define CPPFFT_FFT_N_AS_HPP

#include <iterator>
#include <type_traits>
#include "./detail/is_inputtable_iterator.hpp"
#include "./detail/is_outputtable_iterator.hpp"
#include "./fast_fourier_transform.hpp"

namespace cppfft {

template <
    typename ComplexType,
    typename InputIterator,
    typename DifferenceType,
    typename OutputIterator>
inline auto fft_n_as(InputIterator first, DifferenceType size, OutputIterator result)
    -> std::enable_if_t<
        cppfft::detail::is_inputtable_iterator_v<InputIterator>
            && cppfft::detail::is_outputtable_iterator_v<OutputIterator>,
        OutputIterator>
{
    return cppfft::fast_fourier_transform<ComplexType, DifferenceType>(size)
        .forward(first, result);
}

} // namespace cppfft

#endif // #ifndef CPPFFT_FFT_N_AS_HPP
