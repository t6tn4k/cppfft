#ifndef CPPFFT_FFT_AS_HPP
#define CPPFFT_FFT_AS_HPP

#include <iterator>
#include <type_traits>
#include "./detail/is_forwardable_iterator.hpp"
#include "./detail/is_outputtable_iterator.hpp"
#include "./fast_fourier_transform.hpp"

namespace cppfft {

template <
    typename ComplexType,
    typename ForwardIterator1,
    typename ForwardIterator2,
    typename OutputIterator>
inline auto fft_as(ForwardIterator1 first, ForwardIterator2 last, OutputIterator result)
    -> std::enable_if_t<
        cppfft::detail::is_forwardable_iterator_v<ForwardIterator1>
            && cppfft::detail::is_outputtable_iterator_v<OutputIterator>,
        OutputIterator>
{
    using difference_type = typename std::iterator_traits<ForwardIterator1>::difference_type;
    return cppfft::fast_fourier_transform<ComplexType, difference_type>(std::distance(first, last))
        .forward(first, result);
}

} // namespace cppfft

#endif // #ifndef CPPFFT_FFT_AS_HPP
