cppfft
========

C++ mixed-radix fast Fourier transform

## Synopsis

```cpp
namespace cppfft {

constexpr struct forward_t {} forward{};
constexpr struct inverse_t {} inverse{};

template <typename ComplexType, typename DifferenceType = std::ptrdiff_t>
class fast_fourier_transform
{
public:
    using value_type = ComplexType;
    using element_type = typename value_type::value_type;
    using difference_type = DifferenceType;

    explicit fast_fourier_transform(difference_type n);

    auto size() const noexcept -> difference_type const&;

    template <typename TransformType, typename InputIterator, typename OutputIterator>
    auto operator()(TransformType&&, InputIterator first, OutputIterator result) -> OutputIterator;

    template <typename InputIterator, typename OutputIterator>
    auto operator()(InputIterator first, OutputIterator result) -> OutputIterator;

    template <typename InputIterator, typename OutputIterator>
    auto forward(InputIterator first, OutputIterator result) -> OutputIterator;

    template <typename InputIterator, typename OutputIterator>
    auto inverse(InputIterator first, OutputIterator result) -> OutputIterator;
};

template <typename ForwardIterator1, typename ForwardIterator2, typename OutputIterator>
auto fft(ForwardIterator1 first, ForwardIterator2 last, OutputIterator result) -> OutputIterator;

template <typename ComplexType, typename ForwardIterator1, typename ForwardIterator2, typename OutputIterator>
auto fft_as(ForwardIterator1 first, ForwardIterator2 last, OutputIterator result) -> OutputIterator;

template <typename InputIterator, typename DifferenceType, typename OutputIterator>
auto fft_n(InputIterator first, DifferenceType size, OutputIterator result) -> OutputIterator;

template <typename ComplexType, typename InputIterator, typename DifferenceType, typename OutputIterator>
auto fft_n_as(InputIterator first, DifferenceType size, OutputIterator result) -> OutputIterator;

template <typename ForwardIterator1, typename ForwardIterator2, typename OutputIterator>
auto ifft(ForwardIterator1 first, ForwardIterator2 last, OutputIterator result) -> OutputIterator;

template <typename ComplexType, typename ForwardIterator1, typename ForwardIterator2, typename OutputIterator>
auto ifft_as(ForwardIterator1 first, ForwardIterator2 last, OutputIterator result) -> OutputIterator;

template <typename InputIterator, typename DifferenceType, typename OutputIterator>
auto ifft_n(InputIterator first, DifferenceType size, OutputIterator result) -> OutputIterator;

template <typename ComplexType, typename InputIterator, typename DifferenceType, typename OutputIterator>
auto ifft_n_as(InputIterator first, DifferenceType size, OutputIterator result) -> OutputIterator;

} // namespace cppfft
```
