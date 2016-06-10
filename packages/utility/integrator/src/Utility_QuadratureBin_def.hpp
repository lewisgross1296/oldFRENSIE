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
    d_integral( RQT::zero() ),
    d_integral_abs( RQT::zero() ),
    d_integral_asc( RQT::zero() ),
    d_absolute_error( RQT::zero() )
{ /* ... */ }

// Constructor
template<typename IndepQuantity, typename ResultQuantity>
QuadratureBin<IndepQuantity,ResultQuantity>::QuadratureBin(
                                             const IndepQuantity& lower_limit,
                                             const IndepQuantity& upper_limit )
  : d_lower_limit( lower_limit ),
    d_upper_limit( upper_limit ),
    d_integral( RQT::zero() ),
    d_integral_abs( RQT::zero() ),
    d_integral_asc( RQT::zero() ),
    d_absolute_error( RQT::zero() )
{
  // Make sure the limits are valid
  testPrecondition( lower_limit <= upper_limit );
}

// Copy constructor
template<typename IndepQuantity, typename ResultQuantity>
QuadratureBin<IndepQuantity,ResultQuantity>::QuadratureBin(
                                               const QuadratureBin& other_bin )
  : d_lower_limit( other_bin.d_lower_limit ),
    d_upper_limit( other_bin.d_upper_limit ),
    d_integral( other_bin.d_integral ),
    d_integral_abs( other_bin.d_integral_abs ),
    d_integral_asc( other_bin.d_integral_asc ),
    d_absolute_error( other_bin.d_absolute_error )
{
  // Make sure the limits are valid
  testPrecondition( other_bin.d_lower_limit <= other_bin.d_upper_limit );
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
    d_integral = other_bin.d_integral;
    d_integral_abs = other_bin.d_integral_abs;
    d_integral_asc = other_bin.d_integral_asc;
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
const IndepQuantity&
QuadratureBin<IndepQuantity,ResultQuantity>::getLowerLimit() const
{
  return d_lower_limit;
}

// Get the upper limit of the bin
template<typename IndepQuantity, typename ResultQuantity>
const IndepQuantity&
QuadratureBin<IndepQuantity,ResultQuantity>::getUpperLimit() const
{
  return d_upper_limit;
}

// Set the integral of the bin
template<typename IndepQuantity, typename ResultQuantity>
void QuadratureBin<IndepQuantity,ResultQuantity>::setIntegral(
                                               const ResultQuantity& integral )
{
  // Make sure the result is valid
  testPrecondition( !RQT::isnaninf( integral ) );

  d_integral = integral;
}

// Get the integral of the bin
template<typename IndepQuantity, typename ResultQuantity>
const ResultQuantity&
QuadratureBin<IndepQuantity,ResultQuantity>::getIntegral() const
{
  return d_integral;
}

// Get the integral of the bin
template<typename IndepQuantity, typename ResultQuantity>
ResultQuantity& QuadratureBin<IndepQuantity,ResultQuantity>::getIntegral()
{
  return d_integral;
}

// Set the integral abs of the bin
template<typename IndepQuantity, typename ResultQuantity>
void QuadratureBin<IndepQuantity,ResultQuantity>::setIntegralAbs(
                                           const ResultQuantity& integral_abs )
{
  // Make sure the integral abs is valid
  testPrecondition( !RQT::isnaninf( integral_abs ) );

  d_integral_abs = integral_abs;
}

// Get the integral abs of the bin
template<typename IndepQuantity, typename ResultQuantity>
const ResultQuantity&
QuadratureBin<IndepQuantity,ResultQuantity>::getIntegralAbs() const
{
  return d_integral_abs;
}

// Get the integral abs of the bin
template<typename IndepQuantity, typename ResultQuantity>
ResultQuantity& QuadratureBin<IndepQuantity,ResultQuantity>::getIntegralAbs()
{
  return d_integral_abs;
}

// Set the integral asc of the bin
template<typename IndepQuantity, typename ResultQuantity>
void QuadratureBin<IndepQuantity,ResultQuantity>::setIntegralAsc(
                                           const ResultQuantity& integral_asc )
{
  // Make sure the integral asc is valid
  testPrecondition( !RQT::isnaninf( integral_asc ) );

  d_integral_asc = integral_asc;
}

// Get the integral asc of the bin
template<typename IndepQuantity, typename ResultQuantity>
const ResultQuantity&
QuadratureBin<IndepQuantity,ResultQuantity>::getIntegralAsc() const
{
  return d_integral_asc;
}

// Get the integral asc of the bin
template<typename IndepQuantity, typename ResultQuantity>
ResultQuantity& QuadratureBin<IndepQuantity,ResultQuantity>::getIntegralAsc()
{
  return d_integral_asc;
}

// Set the absolute error of the integral
template<typename IndepQuantity, typename ResultQuantity>
void QuadratureBin<IndepQuantity,ResultQuantity>::setAbsoluteError(
                                                  const ResultQuantity& error )
{
  // Make sure the error is valid
  testPrecondition( !RQT::isnaninf( error ) );

  d_absolute_error = error;
}

// Get the absolute error of the integral
template<typename IndepQuantity, typename ResultQuantity>
const ResultQuantity&
QuadratureBin<IndepQuantity,ResultQuantity>::getAbsoluteError() const
{
  return d_absolute_error;
}

// Get the absolute error of the integral
template<typename IndepQuantity, typename ResultQuantity>
ResultQuantity& QuadratureBin<IndepQuantity,ResultQuantity>::getAbsoluteError()
{
  return d_absolute_error;
}

// Set the result and error of the bin
template<typename IndepQuantity, typename ResultQuantity>
void QuadratureBin<IndepQuantity,ResultQuantity>::setIntegralData(
                                            const ResultQuantity& result,
                                            const ResultQuantity& integral_abs,
                                            const ResultQuantity& integral_asc,
                                            const ResultQuantity error )
{
  // Make sure the integral is valid
  testPrecondition( !RQT::isnaninf( integral ) );
  // Make sure the integral abs is valid
  testPrecondition( !RQT::isnaninf( integral_abs ) );
  // Make sure the integral asc is valid
  testPrecondition( !RQT::isnaninf( integral_asc ) );
  // Make sure the error is valid
  testPrecondition( !RQT::isnaninf( error ) );

  this->setIntegral( integral );
  this->setIntegralAbs( integral_abs );
  this->setIntegralAsc( integral_asc );
  this->setAbsoluteError( error );
}
  
} // end Utility namespace

#endif // end UTILITY_QUADRATURE_BIN_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_QuadratureBin_def.hpp
//---------------------------------------------------------------------------//
