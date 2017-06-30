//---------------------------------------------------------------------------//
//!
//! \file   Utility_ArrayView_def.hpp
//! \author Alex Robinson
//! \brief  The array view class definition
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_ARRAY_VIEW_DEF_HPP
#define UTILITY_ARRAY_VIEW_DEF_HPP

namespace Utility{

// Default constructor
template<typename T>
ArrayView<T>::ArrayView()
  : View<T*>()
{ /* ... */ }

// Iterator constructor
template<typename T>
ArrayView<T>::ArrayView( T* start, T* end )
  : View<T*>( start, end )
{ /* ... */ }

// Range constructor
template<typename T>
ArrayView<T>::ArrayView( T* array_start,
                         const typename ArrayView<T>::size_type array_size )
  : View<T*>( array_start, array_start + array_size )
{ /* ... */ }

// Vector constructor
template<typename T>
ArrayView<T>::ArrayView( std::vector<T>& vector )
  : ArrayView<T>( vector.data(), vector.size() )
{ /* ... */ }

// Const vector constructor
template<typename T>
template<typename U>
ArrayView<T>::ArrayView( const std::vector<U>& vector )
  : ArrayView<T>( vector.data(), vector.size() )
{ /* ... */ }

// Array constructor
template<typename T>
template<size_t N>
ArrayView<T>::ArrayView( std::array<T,N>& array )
  : ArrayView<T>( array.data(), N )
{ /* ... */ }

// Const array constructor
template<typename T>
template<typename U, size_t N>
ArrayView<T>::ArrayView( const std::array<U,N>& array )
  : ArrayView<T>( array.data(), N )
{ /* ... */ }

// Copy constructor
template<typename T>
ArrayView<T>::ArrayView( ArrayView<T>& other_view )
  : View<T*>( other_view )
{ /* ... */ }

// Const array view copy constructor
template<typename T>
template<typename U>
ArrayView<T>::ArrayView( const ArrayView<U>& other_view )
  : View<T*>( other_view )
{ /* ... */ }

// Assignment operator
template<typename T>
ArrayView<T>& ArrayView<T>::operator=( ArrayView<T>& other_view )
{
  if( this != &other_view )
    View<T*>::operator=( other_view );

  return *this;
}

// Const view assignment operator
template<typename T>
template<typename U>
ArrayView<T>& ArrayView<T>::operator=( const ArrayView<U>& other_view )
{
  if( this != &other_view )
    View<T*>::operator=( other_view );

  return *this;
}

// Return a sub-array view
template<typename T>
ArrayView<T> ArrayView<T>::operator()(
                            const typename ArrayView<T>::size_type offset,
                            const typename ArrayView<T>::size_type size ) const
{
  // Make sure that the offset is valid
  testPrecondition( offset < this->size() );
  // Make sure that the slice is valid
  testPrecondition( offset + size <= this->size() );
  
  return ArrayView<T>( const_cast<T*>(this->begin())+offset, size );
}

// Return a sub-array view
template<typename T>
ArrayView<T> ArrayView<T>::operator|( const Utility::Slice& slice ) const
{
  // Make sure that the slice offset is valid
  testPrecondition( slice.offset() < this->size() );
  // Make sure that the slice is valid
  testPrecondition( slice.offset() + slice.extent() <= this->size() );

  return ArrayView<T>( const_cast<T*>(this->begin())+slice.offset(),
                       slice.extent() );
}

// Return a const array view
template<typename T>
ArrayView<const typename std::remove_const<T>::type> ArrayView<T>::toConst() const
{
  return ArrayView<const typename std::remove_const<T>::type>( this->begin(), this->end() );
}

// Implicitly convert to a const array view
template<typename T>
ArrayView<T>::operator ArrayView<const typename std::remove_const<T>::type>() const
{
  return this->toConst();
}
  
} // end Utility namespace

#endif // end UTILITY_ARRAY_VIEW_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_ArrayView_def.hpp
//---------------------------------------------------------------------------//

