#ifndef H_FS_CONCATENATE
#define H_FS_CONCATENATE

#include <type_traits>

namespace flexiblesusy {
namespace detail
{
template<class T> 
struct value_type
{
  using type = typename T::value_type;
};

template<class ...Args>
struct common_value_type
{
  using type = typename std::common_type<
    typename value_type<Args>::type...
  >::type;
};

template<int ...Args> struct sum;

template<int V1, int ...Tail>
struct sum<V1, Tail...>
{
  static constexpr int value = V1 + sum<Tail...>::value;
};

template<> struct sum<>
{
  static constexpr int value = 0;
};

namespace result_of
{
template<class ...Args>
struct concatenate
{
  using type = std::array<
    typename common_value_type<
      typename std::decay<Args>::type...
    >::type,
    sum<std::tuple_size<
      typename std::decay<Args>::type
    >::value...>::value
  >;
};
}

template<class ...Args>
struct concatenate_impl;

template<class T, class ...Args>
struct concatenate_impl<T, Args...>
{
  template<class OutputIterator>
  void operator()( OutputIterator it, T &&t, Args &&...args )
  {
    auto next = std::copy( t.begin(), t.end(), it );
    concatenate_impl<Args...>{}( next, std::forward<Args>( args )... );
  }
};

template<>
struct concatenate_impl<>
{
  template<class OutputIterator>
  void operator()( OutputIterator it )
  {}
};
}

template<class ...Args>
typename detail::result_of::concatenate<Args...>::type
concatenate( Args &&...args )
{
  using result_type =
    typename detail::result_of::concatenate<Args...>::type;
  result_type result;  

  detail::concatenate_impl<Args...>{}( result.begin(),
    std::forward<Args>( args )... );

  return result;
}
}

#endif

