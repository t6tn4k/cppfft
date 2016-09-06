#ifndef CPPFFT_DETAIL_IS_OUTPUTTABLE_ITERATOR_HPP
#define CPPFFT_DETAIL_IS_OUTPUTTABLE_ITERATOR_HPP

#include <iterator>
#include <type_traits>

namespace cppfft { namespace detail {

template <typename T>
struct is_outputtable_iterator
    : std::integral_constant<bool,
        std::is_same<
            typename std::iterator_traits<T>::iterator_category,
            std::output_iterator_tag>::value
            || std::is_convertible<
                typename std::iterator_traits<T>::iterator_category,
                std::forward_iterator_tag>::value>
{
};

template <typename T>
constexpr auto is_outputtable_iterator_v = cppfft::detail::is_outputtable_iterator<T>::value;

} } // namespace cppfft::detail

#endif // #ifndef CPPFFT_DETAIL_IS_OUTPUTTABLE_ITERATOR_HPP
