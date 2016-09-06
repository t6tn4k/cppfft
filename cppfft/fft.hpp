#ifndef CPPFFT_FFT_HPP
#define CPPFFT_FFT_HPP

#include <iterator>
#include <type_traits>
#include "./detail/is_forwardable_iterator.hpp"
#include "./fast_fourier_transform.hpp"

namespace cppfft {

template <typename ForwardIterator1, typename ForwardIterator2, typename OutputIterator>
inline auto fft(ForwardIterator1 first, ForwardIterator2 last, OutputIterator result)
    -> std::enable_if_t<
        cppfft::detail::is_forwardable_iterator_v<ForwardIterator1>,
        OutputIterator>
{
    using value_type = typename std::iterator_traits<ForwardIterator1>::value_type;
    using difference_type = typename std::iterator_traits<ForwardIterator1>::difference_type;
    return cppfft::fast_fourier_transform<value_type, difference_type>(std::distance(first, last))
        .forward(first, result);
}

} // namespace cppfft

#endif // #ifndef CPPFFT_FFT_HPP
