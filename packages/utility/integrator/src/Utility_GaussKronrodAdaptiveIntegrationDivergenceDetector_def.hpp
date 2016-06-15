//---------------------------------------------------------------------------//
//!
//! \file   Utility_GaussKronrodAdaptiveIntegrationDivergenceDetector_def.hpp
//! \author Alex Robinson
//! \brief  The Gauss-Kronrod adaptive int. divergence detector template defs
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_GAUSS_KRONROD_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_DEF_HPP
#define UTILITY_GAUSS_KRONROD_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_DEF_HPP

// FRENSIE Includes
#include "Utility_IntegratorException.hpp"
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_ContractException.hpp"

namespace Utility{

// Constructor
template<typename IntegralQuantity>
GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::GaussKronrodAdaptiveIntegrationDivergenceDetector(
                    const unsigned slowly_converging_integral_count_limit,
                    const unsigned diverging_integral_count_limit,
                    const typename IQT::RawType integral_change_upper_limit,
                    const typename IQT::RawType error_slow_change_lower_limit,
                    const unsigned iteration_limit )
  : AdaptiveIntegralDivergenceDetector<IntegralQuantity(
                                                 integral_change_upper_limit,
                                                 error_slow_change_lower_limit,
                                                 interation_limit ),
    d_slowly_converging_integral_counter( 0u ),
    d_diverging_integral_counter( 0u ),
    d_slowly_converging_integral_count_limit(
                                      slowly_converging_integral_count_limit ),
    d_diverging_integral_count_limit( diverging_integral_count_limit )
{
  // Make sure the slowly converging integral count limit is valid
  testPrecondition( slowly_converging_integral_count_limit > 0u );
  // Make sure the diverging integral count limit is valid
  testPrecondition( diverging_integral_count_limit > 0u );
}

// Check for a diverging integral estimate
/*! \details If divergence is detected an IntegratorException will be thrown.
 */
template<typename IntegralQuantity>
void GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::checkForDivergence(
                                     const IntegralQuantity& original_integral,
                                     const IntegralQuantity& original_error,
                                     const IntegralQuantity& refined_integral,
                                     const IntegralQuantity& refined_error,
                                     const unsigned number_of_iterations )
{
  // Check if the slowly converging integral counter needs to be incremented
  if( this->isIntegralVerySlowlyConverging( original_integral,
                                            original_error,
                                            refined_integral,
                                            refined_error ) )
  {
    this->incrementSlowlyConvergingIntegralCounter();
  }

  // Check if the diverging integral counter needs to be incremented
  if( this->isIntegralDiverging( original_error,
                                 refined_error,
                                 number_of_iterations ) )
  {
    this->incrementDivergingIntegralCounter();
  }

  // Check if the slowly converging integral count limit has been reached
  if( this->hasSlowlyConvergingIntegralCountLimitBeenReached() )
  {
    THROW_EXCEPTION( Utility::IntegratorException,
                     "Error: Roundoff error prevented the integral from "
                     "converging (" << number_of_iterations <<
                     " iterations)!" );
  }

  // Check if the diverging integral count limit has been reached
  if( this->hasDivergingIntegralCountLimitBeenReached() )
  {
    THROW_EXCEPTION( Utility::IntegratorException,
                     "Error: Roundoff error prevented the integral from "
                     "converging (" << number of iterations <<
                     " iterations)!" );
  }
}

// Get the slowly converging integral count
template<typename IntegralQuantity>
unsigned GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::getSlowlyConvergingIntegralCount() const
{
  return d_slowly_converging_integral_counter;
}

// Get the diverging integral count
template<typename IntegralQuantity>
unsigned GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::getDivergingIntegralCount() const
{
  return d_diverging_integral_counter;
}

// Get the slowly converging integral count limit
template<typename IntegralQuantity>
unsigned GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::getSlowlyConvergingIntegralCountLimit() const
{
  return d_slowly_converging_integral_count_limit;
}

// Get the diverging integral count limit
template<typename IntegralQuantity>
unsigned GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::getDivergingIntegralCountLimit() const
{
  return d_diverging_integral_count_limit;
}

// Increment the slowly converging integral counter
template<typename IntegralQuantity>
void GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::incrementSlowlyConvergingIntegralCounter()
{
  ++d_slowly_converging_integral_counter;
}

// Increment the diverging integral counter
template<typename IntegralQuantity>
void GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::incrementDivergingIntegralCounter()
{
  ++d_diverging_integral_counter;
}

// Test if the slowly converging integral count limit has been reached
template<typename IntegralQuantity>
bool GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::hasSlowlyConvergingIntegralCountLimitBeenReached() const
{
  return d_slowly_converging_integral_counter >=
    d_slowly_converging_integral_count_limit;
}

// Test if the diverging integral count limit has been reached
template<typename IntegralQuantity>
bool GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::hasDivergingIntegralCountLimitBeenReached() const
{
  return d_diverging_integral_counter >=
    d_diverging_integral_count_limit;
}
  
} // end Utility namespace

#endif // end UTILITY_GAUSS_KRONROD_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_GaussKronrodAdaptiveIntegrationDivergenceDetector_def.hpp
//---------------------------------------------------------------------------//
