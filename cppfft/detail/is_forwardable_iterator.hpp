#ifndef CPPFFT_DETAIL_IS_FORWARDABLE_ITERATOR_HPP
#define CPPFFT_DETAIL_IS_FORWARDABLE_ITERATOR_HPP

#include <iterator>
#include <type_traits>

namespace cppfft { namespace detail {

template <typename T>
struct is_forwardable_iterator
    : std::is_convertible<
        typename std::iterator_traits<T>::iterator_category,
        std::forward_iterator_tag>
{
};

template <typename T>
constexpr auto is_forwardable_iterator_v = cppfft::detail::is_forwardable_iterator<T>::value;

} } // namespace cppfft::detail

#endif // #ifndef CPPFFT_DETAIL_IS_FORWARDABLE_ITERATOR_HPP
