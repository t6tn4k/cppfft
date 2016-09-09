#ifndef CPPFFT_IFFT_N_HPP
#define CPPFFT_IFFT_N_HPP

#include <iterator>
#include <type_traits>
#include "./detail/is_inputtable_iterator.hpp"
#include "./detail/is_outputtable_iterator.hpp"
#include "./fast_fourier_transform.hpp"

namespace cppfft {

template <typename InputIterator, typename DifferenceType, typename OutputIterator>
inline auto ifft_n(InputIterator first, DifferenceType size, OutputIterator result)
    -> std::enable_if_t<
        cppfft::detail::is_inputtable_iterator_v<InputIterator>
            && cppfft::detail::is_outputtable_iterator_v<OutputIterator>,
        OutputIterator>
{
    using value_type = typename std::iterator_traits<InputIterator>::value_type;
    return cppfft::fast_fourier_transform<value_type, DifferenceType>(size)
        .inverse(first, result);
}

} // namespace cppfft

#endif // #ifndef CPPFFT_IFFT_N_HPP
