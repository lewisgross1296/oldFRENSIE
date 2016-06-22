//---------------------------------------------------------------------------//
//!
//! \file   Utility_WynnEpsilonPriorityQueue.hpp
//! \author Alex Robinson
//! \brief  The Wynn-Epsilon priority queue (heap) declaration
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_WYNN_EPSILON_PRIORITY_QUEUE_HPP
#define UTILITY_WYNN_EPSILON_PRIORITY_QUEUE_HPP

// Std Lib Includes
#include <vector>
#include <function>

// FRENSIE Includes
#include "QuadratureBin.hpp"

namespace Utility{

//! The Wynn-Epsilon priority queue (heap) container adaptor
template<typename ArgQuantity,
         typename IntegralQuantity,
         typename Sequence = std::vector<QuadratureBin<ArgQuantity,IntegralQuantity> > >
class WynnEpsilonPriorityQueue
{

public:

  //! The value type (std lib compliant)
  typedef typename Sequence::value_type value_type;

  //! The reference type (std lib compliant)
  typedef typename Sequence::reference reference;

  //! The const reference type (std lib compliant)
  typedef typename Sequence::const_reference const_reference;

  //! The size type (std lib compliant)
  typedef typename Sequence::size_type size_type

  //! The container type (std lib compliant)
  typedef Sequence container_type;

  //! Constructor
  explicit WynnEpsilonPriorityQueue( const size_type max_bisection_level,
                                     const Sequence& sequence = Sequence() );

  //! Constructor
  explicit WynnEpsilonPriorityQueue( const size_type max_bisection_level,
                                     const Sequence&& sequence = Sequence() );

  //! Constructor
  template<typename InputIterator>
  explicit WynnEpsilonPriorityQueue( const size_type max_bisection_level,
                                     InputIterator first,
                                     InputIterator last,
                                     const Sequence& sequence = Sequence() );

  //! Constructor
  template<typename InputIterator>
  explicit WynnEpsilonPriorityQueue( const size_type max_bisection_level,
                                     InputIterator first,
                                     InputIterator last,
                                     const Sequence&& sequence = Sequence() );

  //! Test if the queue is empty
  bool empty() const;

  //! Get the number of elements in the queue
  size_type size() const;

  //! Get the first element in the queue (read-only)
  const_reference top() const;

  //! Insert an element into the queue
  void push( const value_type& new_element );

  //! Insert an element into the queue
  void push( const value_type&& new_element );

  //! Construct and insert new elements into the queue
  template<typename... ElementArgs>
  void emplace( ElementArgs&&... element_args );

  //! Remove the first element from the queue
  void pop();

  //! Set the max bisection level for bins
  void setMaxBisectionLevel( const size_type bisection_level );

  //! Increment the max bisection level for bins (by one)
  void incrementBisectionLevel();

  //! Get the max bisection level for bins
  size_type getMaxBisectionLevel() const;

private:

  // The max bisection level
  size_type d_max_bisection_level;

  // The sequence
  Sequence d_sequence;

  // The comparison function
  std::function<bool(const_reference,const_reference)> d_compare_fn;
};
  
} // end Utility namespace

//---------------------------------------------------------------------------//
// Template Includes
//---------------------------------------------------------------------------//

#include "Utility_WynnEpsilonPriorityQueue_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_WYNN_EPSILON_PRIORITY_QUEUE_HPP

//---------------------------------------------------------------------------//
// end Utility_WynnEpsilonPriorityQueue.hpp
//---------------------------------------------------------------------------//
