#ifndef CPPFFT_IFFT_AS_HPP
#define CPPFFT_IFFT_AS_HPP

#include <iterator>
#include <type_traits>
#include "./detail/is_forwardable_iterator.hpp"
#include "./fast_fourier_transform.hpp"

namespace cppfft {

template <
    typename ComplexType,
    typename ForwardIterator1,
    typename ForwardIterator2,
    typename OutputIterator>
inline auto ifft_as(ForwardIterator1 first, ForwardIterator2 last, OutputIterator result)
    -> std::enable_if_t<
        cppfft::detail::is_forwardable_iterator_v<ForwardIterator1>,
        OutputIterator>
{
    using difference_type = typename std::iterator_traits<ForwardIterator1>::difference_type;
    return cppfft::fast_fourier_transform<ComplexType, difference_type>(std::distance(first, last))
        .inverse(first, result);
}

} // namespace cppfft

#endif // #ifndef CPPFFT_IFFT_AS_HPP
