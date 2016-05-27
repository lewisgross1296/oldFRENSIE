//---------------------------------------------------------------------------//
//!
//! \file   Utility_QuadratureBin_def.hpp
//! \author Alex Robinson
//! \brief  The quadrature bin definition
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_QUADRATURE_BIN_DEF_HPP
#define UTILITY_QUADRATURE_BIN_DEF_HPP

// FRENSIE Includes
#include "Utility_ContractException.hpp"

namespace Utility{

// Default constructor
template<typename IndepQuantity, typename ResultQuantity>
QuadratureBin<IndepQuantity,ResultQuantity>::QuadratureBin()
  : d_lower_limit( IQT::zero() ),
    d_upper_limit( IQT::zero() ),
    d_result( RQT::zero() ),
    d_absolute_error( RQT::zero() )
{ /* ... */ }

// Constructor
template<typename IndepQuantity, typename ResultQuantity>
QuadratureBin<IndepQuantity,ResultQuantity>::QuadratureBin(
                                              const IndepQuantity lower_limit,
                                              const IndepQuantity upper_limit )
  : d_lower_limit( lower_limit ),
    d_upper_limit( upper_limit ),
    d_result( RQT::zero() ),
    d_absolute_error( RQT::zero() )
{
  // Make sure the limits are valid
  testPrecondition( lower_limit <= upper_limit );
}

// Assignment operator
template<typename IndepQuantity, typename ResultQuantity>
QuadratureBin<IndepQuantity,ResultQuantity>&
QuadratureBin<IndepQuantity,ResultQuantity>::operator=(
                 const QuadratureBin<IndepQuantity,ResultQuantity>& other_bin )
{
  if( this != &other_bin )
  {
    d_lower_limit = other_bin.d_lower_limit;
    d_upper_limit = other_bin.d_upper_limit;
    d_result = other_bin.d_result;
    d_absolute_error = other_bin.d_absolute_error;
  }

  return *this;
}

// Comparison operator
/*! \details The absolute errors will be used to do the comparison
 */
template<typename IndepQuantity, typename ResultQuantity>
bool QuadratureBin<IndepQuantity,ResultQuantity>::operator<(
                 const QuadratureBin<IndepQuantity,ResultQuantity>& other_bin )
{
  return d_absolute_error < other_bin.absolute_error;
}

// Get the lower limit of the bin
template<typename IndepQuantity, typename ResultQuantity>
IndepQuantity
QuadratureBin<IndepQuantity,ResultQuantity>::getLowerLimit() const
{
  return d_lower_limit;
}

// Get the upper limit of the bin
template<typename IndepQuantity, typename ResultQuantity>
IndepQuantity
QuadratureBin<IndepQuantity,ResultQuantity>::getUpperLimit() const
{
  return d_upper_limit;
}

// Set the result of the bin
template<typename IndepQuantity, typename ResultQuantity>
void QuadratureBin<IndepQuantity,ResultQuantity>::setResult(
                                                 const ResultQuantity result )
{
  // Make sure the result is valid
  testPrecondition( !RQT::isnaninf( result ) );

  d_result = result;
}

// Get the result of the bin
template<typename IndepQuantity, typename ResultQuantity>
ResultQuantity QuadratureBin<IndepQuantity,ResultQuantity>::getResult() const
{
  return d_result;
}

// Set the error of the result
template<typename IndepQuantity, typename ResultQuantity>
void QuadratureBin<IndepQuantity,ResultQuantity>::setAbsoluteError(
                                                   const ResultQuantity error )
{
  // Make sure the error is valid
  testPrecondition( !RQT::isnaninf( error ) );

  d_absolute_error = error;
}

// Get the error of the result
template<typename IndepQuantity, typename ResultQuantity>
ResultQuantity
QuadratureBin<IndepQuantity,ResultQuantity>::getAbsoluteError() const
{

}

// Set the result and error of the bin
template<typename IndepQuantity, typename ResultQuantity>
void QuadratureBin<IndepQuantity,ResultQuantity>::setResultAndAbsoluteError(
                                                   const ResultQuantity result,
                                                   const ResultQuantity error )
{
  // Make sure the result is valid
  testPrecondition( !RQT::isnaninf( result ) );
  // Make sure the error is valid
  testPrecondition( !RQT::isnaninf( error ) );

  this->setResult( result );
  this->setAbsoluteError( error );
}
  
} // end Utility namespace

#endif // end UTILITY_QUADRATURE_BIN_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_QuadratureBin_def.hpp
//---------------------------------------------------------------------------//
