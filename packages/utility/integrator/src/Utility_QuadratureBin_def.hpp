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
template<typename ArgQuantity, typename IntegralQuantity>
QuadratureBin<ArgQuantity,IntegralQuantity>::QuadratureBin()
  : d_lower_limit( AQT::zero() ),
    d_upper_limit( AQT::zero() ),
    d_integral( IQT::zero() ),
    d_integral_abs( IQT::zero() ),
    d_integral_asc( IQT::zero() ),
    d_absolute_error( IQT::zero() ),
    d_bisection_level( 0u )
{ /* ... */ }

// Constructor
template<typename ArgQuantity, typename IntegralQuantity>
QuadratureBin<ArgQuantity,IntegralQuantity>::QuadratureBin(
                                              const ArgQuantity& lower_limit,
                                              const ArgQuantity& upper_limit )
  : d_lower_limit( lower_limit ),
    d_upper_limit( upper_limit ),
    d_integral( IQT::zero() ),
    d_integral_abs( IQT::zero() ),
    d_integral_asc( IQT::zero() ),
    d_absolute_error( IQT::zero() ),
    d_bisection_level( 0u )
{
  // Make sure the limits are valid
  testPrecondition( lower_limit <= upper_limit );
}

// Copy constructor
template<typename ArgQuantity, typename IntegralQuantity>
QuadratureBin<ArgQuantity,IntegralQuantity>::QuadratureBin(
                                               const QuadratureBin& other_bin )
  : d_lower_limit( other_bin.d_lower_limit ),
    d_upper_limit( other_bin.d_upper_limit ),
    d_integral( other_bin.d_integral ),
    d_integral_abs( other_bin.d_integral_abs ),
    d_integral_asc( other_bin.d_integral_asc ),
    d_absolute_error( other_bin.d_absolute_error ),
    d_bisection_level( other_bin.d_bisection_limit )
{
  // Make sure the limits are valid
  testPrecondition( other_bin.d_lower_limit <= other_bin.d_upper_limit );
}

// Assignment operator
template<typename ArgQuantity, typename IntegralQuantity>
QuadratureBin<ArgQuantity,IntegralQuantity>&
QuadratureBin<ArgQuantity,IntegralQuantity>::operator=(
                 const QuadratureBin<ArgQuantity,IntegralQuantity>& other_bin )
{
  if( this != &other_bin )
  {
    d_lower_limit = other_bin.d_lower_limit;
    d_upper_limit = other_bin.d_upper_limit;
    d_integral = other_bin.d_integral;
    d_integral_abs = other_bin.d_integral_abs;
    d_integral_asc = other_bin.d_integral_asc;
    d_absolute_error = other_bin.d_absolute_error;
    d_bisection_limit = other_bin.d_bisection_limit;
  }

  return *this;
}

// Get the lower limit of the bin
template<typename ArgQuantity, typename IntegralQuantity>
const ArgQuantity&
QuadratureBin<ArgQuantity,IntegralQuantity>::getLowerLimit() const
{
  return d_lower_limit;
}

// Get the upper limit of the bin
template<typename ArgQuantity, typename IntegralQuantity>
const ArgQuantity&
QuadratureBin<ArgQuantity,IntegralQuantity>::getUpperLimit() const
{
  return d_upper_limit;
}

// Set the integral of the bin
template<typename ArgQuantity, typename IntegralQuantity>
void QuadratureBin<ArgQuantity,IntegralQuantity>::setIntegral(
                                             const IntegralQuantity& integral )
{
  // Make sure the result is valid
  testPrecondition( !IQT::isnaninf( integral ) );

  d_integral = integral;
}

// Get the integral of the bin
template<typename ArgQuantity, typename IntegralQuantity>
const IntegralQuantity&
QuadratureBin<ArgQuantity,IntegralQuantity>::getIntegral() const
{
  return d_integral;
}

// Get the integral of the bin
template<typename ArgQuantity, typename IntegralQuantity>
IntegralQuantity& QuadratureBin<ArgQuantity,IntegralQuantity>::getIntegral()
{
  return d_integral;
}

// Set the integral abs of the bin
template<typename ArgQuantity, typename IntegralQuantity>
void QuadratureBin<ArgQuantity,IntegralQuantity>::setIntegralAbs(
                                         const IntegralQuantity& integral_abs )
{
  // Make sure the integral abs is valid
  testPrecondition( !IQT::isnaninf( integral_abs ) );

  d_integral_abs = integral_abs;
}

// Get the integral abs of the bin
template<typename ArgQuantity, typename IntegralQuantity>
const IntegralQuantity&
QuadratureBin<ArgQuantity,IntegralQuantity>::getIntegralAbs() const
{
  return d_integral_abs;
}

// Get the integral abs of the bin
template<typename ArgQuantity, typename IntegralQuantity>
IntegralQuantity& QuadratureBin<ArgQuantity,IntegralQuantity>::getIntegralAbs()
{
  return d_integral_abs;
}

// Set the integral asc of the bin
template<typename ArgQuantity, typename IntegralQuantity>
void QuadratureBin<ArgQuantity,IntegralQuantity>::setIntegralAsc(
                                         const IntegralQuantity& integral_asc )
{
  // Make sure the integral asc is valid
  testPrecondition( !IQT::isnaninf( integral_asc ) );

  d_integral_asc = integral_asc;
}

// Get the integral asc of the bin
template<typename ArgQuantity, typename IntegralQuantity>
const IntegralQuantity&
QuadratureBin<ArgQuantity,IntegralQuantity>::getIntegralAsc() const
{
  return d_integral_asc;
}

// Get the integral asc of the bin
template<typename ArgQuantity, typename IntegralQuantity>
IntegralQuantity& QuadratureBin<ArgQuantity,IntegralQuantity>::getIntegralAsc()
{
  return d_integral_asc;
}

// Set the absolute error of the integral
template<typename ArgQuantity, typename IntegralQuantity>
void QuadratureBin<ArgQuantity,IntegralQuantity>::setAbsoluteError(
                                                const IntegralQuantity& error )
{
  // Make sure the error is valid
  testPrecondition( !IQT::isnaninf( error ) );

  d_absolute_error = error;
}

// Get the absolute error of the integral
template<typename ArgQuantity, typename IntegralQuantity>
const IntegralQuantity&
QuadratureBin<ArgQuantity,IntegralQuantity>::getAbsoluteError() const
{
  return d_absolute_error;
}

// Get the absolute error of the integral
template<typename ArgQuantity, typename IntegralQuantity>
IntegralQuantity&
QuadratureBin<ArgQuantity,IntegralQuantity>::getAbsoluteError()
{
  return d_absolute_error;
}

// Set the result and error of the bin
template<typename ArgQuantity, typename IntegralQuantity>
void QuadratureBin<ArgQuantity,IntegralQuantity>::setIntegralData(
                                          const IntegralQuantity& result,
                                          const IntegralQuantity& integral_abs,
                                          const IntegralQuantity& integral_asc,
                                          const IntegralQuantity error )
{
  // Make sure the integral is valid
  testPrecondition( !IQT::isnaninf( integral ) );
  // Make sure the integral abs is valid
  testPrecondition( !IQT::isnaninf( integral_abs ) );
  // Make sure the integral asc is valid
  testPrecondition( !IQT::isnaninf( integral_asc ) );
  // Make sure the error is valid
  testPrecondition( !IQT::isnaninf( error ) );

  this->setIntegral( integral );
  this->setIntegralAbs( integral_abs );
  this->setIntegralAsc( integral_asc );
  this->setAbsoluteError( error );
}

// Bisect the bin
/*! \details The new bins will have their bisection level increased by one.
 * Any integral data stored in the new bins will also be reset.
 */
template<typename ArgQuantity, typename IntegralQuantity>
void QuadratureBin<ArgQuantity,IntegralQuantity>::bisect(
            QuadratureBin<ArgQuantity,IntegralQuantity>& left_half_bin,
            QuadratureBin<ArgQuantity,IntegralQuantity>& right_half_bin ) const
{
  const ArgQuantity midpoint = (d_lower_limit+d_upper_limit)/2;
  
  left_half_bin.setLimits( d_lower_limit, midpoint );
  left_half_bin.incrementBisectionLevel();
  left_half_bin.resetIntegralData();
  
  right_half_bin.setLimits( midpoint, d_upper_limit );
  right_half_bin.incrementBisectionLevel();
  right_half_bin.resetIntegralData();
}

// Get the bisection level of the bin
/*! \details The bisection level is the number of bisections that were done
 * to get to this bins limits. The first bin created will have a bisection
 * level of 0.
 */
template<typename ArgQuantity, typename IntegralQuantity>
unsigned QuadratureBin<ArgQuantity,IntegralQuantity>::getBisectionLevel() const
{
  return d_bisection_level;
}

// Set the bin limits
template<typename ArgQuantity, typename IntegralQuantity>
void QuadratureBin<ArgQuantity,IntegralQuantity>::setLimits(
                                               const ArgQuantity& lower_limit,
                                               const ArgQuantity& upper_limit )
{
  // Make sure the limits are valid
  testPrecondition( lower_limit <= upper_limit );

  d_lower_limit = lower_limit;
  d_upper_limit = upper_limit;
}

// Increment the bisection level
template<typename ArgQuantity, typename IntegralQuantity>
void QuadratureBin<ArgQuantity,IntegralQuantity>::incrementBisectionLevel()
{
  ++d_bisection_level;
}

// Reset the integral and error of the bin
template<typename ArgQuantity, typename IntegralQuantity>
void QuadratureBin<ArgQuantity,IntegralQuantity>::resetIntegralData()
{
  d_integral = IQT::zero();
  d_integral_abs = IQT::zero();
  d_integral_asc = IQT::zero();
  d_absolute_error = IQT::zero();
}
  
} // end Utility namespace

#endif // end UTILITY_QUADRATURE_BIN_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_QuadratureBin_def.hpp
//---------------------------------------------------------------------------//
