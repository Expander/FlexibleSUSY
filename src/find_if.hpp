#ifndef H_FS_FIND_IF
#define H_FS_FIND_IF

#include <type_traits>
#include <utility>

#include <boost/mpl/next.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/begin_end.hpp>

namespace flexiblesusy {
namespace meta {
namespace detail
{
template<
  class Iterator,
  class End,
  template<typename> class Predicate,
  template<typename> class F,
  class State
> bool find_if_impl( State &&state, std::true_type )
{ return false; }

template<
  class Iterator,
  class End,
  template<typename> class Predicate,
  template<typename> class F,
  class State
> bool find_if_impl( State &&state, std::false_type )
{
  if( Predicate<typename boost::mpl::deref<Iterator>::type>{}(
        state ) )
  {
    F<typename boost::mpl::deref<Iterator>::type>{}( state );
    return true;
  }

  using NextIterator = typename boost::mpl::next<Iterator>::type;

  return find_if_impl<
    NextIterator, End, Predicate, F
  >( std::forward<State>( state ),
     typename std::is_same<NextIterator, End>::type{} );
}
}

template<
  class Sequence,
  template<typename> class Predicate,
  template<typename> class F,
  class State
> bool find_if( State &&state )
{
  return detail::find_if_impl<
    typename boost::mpl::begin<Sequence>::type,
    typename boost::mpl::end<Sequence>::type,
    Predicate, F
  >( std::forward<State>( state ),
     typename std::is_same<
       typename boost::mpl::begin<Sequence>::type,
       typename boost::mpl::end<Sequence>::type
     >::type{} );
}
}
}

#endif

