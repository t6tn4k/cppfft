#ifndef CPPFFT_DETAIL_IS_RANDOM_ACCESSIBLE_ITERATOR_HPP
#define CPPFFT_DETAIL_IS_RANDOM_ACCESSIBLE_ITERATOR_HPP

#include <iterator>
#include <type_traits>

namespace cppfft { namespace detail {

template <typename T>
struct is_random_accessible_iterator
    : std::is_same<
        typename std::iterator_traits<T>::iterator_category,
        std::random_access_iterator_tag>
{
};

template <typename T>
constexpr auto is_random_accessible_iterator_v
    = cppfft::detail::is_random_accessible_iterator<T>::value;

} } // namespace cppfft::detail

#endif // #ifndef CPPFFT_DETAIL_IS_RANDOM_ACCESSIBLE_ITERATOR_HPP
