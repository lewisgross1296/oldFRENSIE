//---------------------------------------------------------------------------//
//!
//! \file   Utility_WynnEpsilonPriorityQueue.hpp
//! \author Alex Robinson
//! \brief  The Wynn-Epsilon priority queue (heap) definition
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_WYNN_EPSILON_PRIORITY_QUEUE_DEF_HPP
#define UTILITY_WYNN_EPSILON_PRIORITY_QUEUE_DEF_HPP

// Std Lib Includes
#include <functional>
#include <algorithm>

namespace Utility{

// Constructor
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::WynnEpsilonPriorityQueue(
                                           const size_type max_bisection_level,
                                           const Sequence& sequence )
  : d_max_bisection_level( max_bisection_level ),
    d_sequence( sequence ),
    d_compare_fn( std::bind<bool>( value_type::compareBisectionLevelsAndBinErrors( _1, _2, max_bisection_level ) ) )
{
  std::make_heap( d_sequence.begin(), d_sequence.end(), d_compare_fn );
}

// Constructor
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::WynnEpsilonPriorityQueue(
                                           const size_type max_bisection_level,
                                           const Sequence&& sequence )
  : d_max_bisection_level( max_bisection_level ),
    d_sequence( std::move( sequence ) ),
    d_compare_fn( std::bind<bool>( value_type::compareBisectionLevelsAndBinErrors( _1, _2, max_bisection_level ) ) )
{
  std::make_heap( d_sequence.begin(), d_sequence.end(), d_compare_fn );
}

// Constructor
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
template<typename InputIterator>
WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::WynnEpsilonPriorityQueue(
                                           const size_type max_bisection_level,
                                           InputIterator first,
                                           InputIterator last,
                                           const Sequence& sequence )
  : d_max_bisection_level( max_bisection_level ),
    d_sequence( sequence ),
    d_compare_fn( std::bind<bool>( value_type::compareBisectionLevelsAndBinErrors( _1, _2, max_bisection_level ) ) )
{
  d_sequence.insert( d_sequence.end(), first, last );
  std::make_heap( d_sequence.begin(), d_sequence.end(), d_compare_fn );
}

// Constructor
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
template<typename InputIterator>
WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::WynnEpsilonPriorityQueue(
                                           const size_type max_bisection_level,
                                           InputIterator first,
                                           InputIterator last,
                                           const Sequence&& sequence )
  : d_max_bisection_level( max_bisection_level ),
    d_sequence( std::move( sequence ) ),
    d_compare_fn( std::bind<bool>( value_type::compareBisectionLevelsAndBinErrors( _1, _2, max_bisection_level ) ) )
{
  d_sequence.insert( d_sequence.end(), first, last );
  std::make_heap( d_sequence.begin(), d_sequence.end(), d_compare_fn );
}

// Test if the queue is empty
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
bool WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::empty() const
{
  return d_sequence.empty();
}

// Get the number of elements in the queue
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
size_type WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::size() const
{
  return d_sequence.size();
}

// Get the first element in the queue (read-only)
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
const_reference WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::top() const
{
  return d_sequence.front();
}

// Insert an element into the queue
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
void WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::push(
                                                const value_type& new_element )
{
  d_sequence.push_back( new_element );
  std::push_heap( d_sequence.begin(), d_sequence.end(), d_compare_fn );
}

// Insert an element into the queue
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
void WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::push(
                                               const value_type&& new_element )
{
  d_sequence.push_back( std::move( new_element ) );
  std::push_heap( d_sequence.begin(), d_sequence.end(), d_compare_fn );
}

// Construct and insert new elements into the queue
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
template<typename... ElementArgs>
void WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::emplace( ElementArgs&&... element_args )
{
  d_sequence.emplace_back( std::forward<ElementArgs>( element_args )... );
  std::push_heap( d_sequence.begin(), d_sequence.end(), d_compare_fn );
}

// Remove the first element from the queue
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
void WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::pop()
{
  std::pop_heap( d_sequence.begin(), d_sequence.end(), d_compare_fn );
  d_sequence.pop_back();
}

// Set the max bisection level for bins
/*! \details The priority queue will be reset after the max bisection level
 * is changed.
 */
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
void WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::setMaxBisectionLevel( const size_type bisection_level )
{
  if( d_max_bisection_level != bisection_level )
  {
    d_max_bisection_level = bisection_level;

    d_compare_fn = std::bind<bool>( value_type::compareBisectionLevelsAndBinErrors( _1, _2, bisection_level ) );
      
    std::make_heap( d_sequence.begin(), d_sequence.end(), d_compare_fn );
  }
}

// Set the max bisection level for bins
/*! \details The priority queue will be reset after the max bisection level
 * is incremented.
 */
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
void WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::incrementMaxBisectionLevel()
{
  ++d_max_bisection_level;

  d_compare_fn = std::bind<bool>( value_type::compareBisectionLevelsAndBinErrors( _1, _2, d_max_bisection_level ) );
      
  std::make_heap( d_sequence.begin(), d_sequence.end(), d_compare_fn );
}

// Get the max bisection level for bins
template<typename ArgQuantity,typename IntegralQuantity,typename Sequence>
size_type WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity,Sequence>::getMaxBisectionLevel() const
{
  return d_max_bisection_level;
}
  
} // end Utility namespace

#endif // end UTILITY_WYNN_EPSILON_PRIORITY_QUEUE_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_WynnEpsilonPriorityQueue_def.hpp
//---------------------------------------------------------------------------//
