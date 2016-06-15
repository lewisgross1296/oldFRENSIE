//---------------------------------------------------------------------------//
//!
//! \file   Utility_AdaptiveIntegrationDivergenceDetector_def.hpp
//! \author Alex Robinson
//! \brief  The adaptive integration divergence detector template defs.
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_DEF_HPP
#define UTILITY_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_DEF_HPP

// FRENSIE Includes
#include "Utility_ContractException.hpp"

namespace Utility{

// Constructor
template<typename IntegralQuantity>
AdaptiveIntegrationDivergenceDetector<IntegralQuantity>::AdaptiveIntegrationDivergenceDetector(
                    const typename IQT::RawType& integral_change_upper_limit,
                    const typename IQT::RawType& error_slow_change_lower_limit,
                    const unsigned iteration_limit )
  : d_integral_change_upper_limit( integral_change_upper_limit ),
    d_error_slow_change_lower_limit( error_slow_change_lower_limit ),
    d_iteration_limit( iteration_limit );
{
  // Make sure the integral change upper limit is valid
  testPrecondition( integral_change_upper_limit < 1.0 );
  testPrecondition( integral_change_upper_limit > 0.0 );
  // Make sure the error slow change lower limit is valid
  testPrecondition( error_slow_change_lower_limit < 1.0 );
  testPrecondition( error_slow_change_lower_limit > 0.0 );
}

// Check for a diverging integral estimator
template<typename IntegralQuantity>
template<typename ArgQuantity>
inline void
AdaptiveIntegrationDivergenceDetector<IntegralQuantity>::checkForDivergence(
             const QuadratureBin<ArgQuantity,IntegralQuantity>& original_bin,
             const QuadratureBin<ArgQuantity,IntegralQuantity>& left_half_bin,
             const QuadratureBin<ArgQuantity,IntegralQuantity>& right_half_bin,
             const unsigned number_of_iterations )
{
  // Calculate the refined integral
  const IntegralQuantity refined_integral =
    left_half_bin.getIntegral()+right_half_bin.getIntegral();

  // Calculate the refined error
  const IntegralQuantity refined_error =
    left_half_bin.getAbsoluteError()+right_half_bin.getAbsoluteError();
  
  this->checkForDivergence( original_bin.getIntegral(),
                            original_bin.getAbsoluteError(),
                            returned_integral,
                            refined_error,
                            number_of_iterations );
}

// Check for a very slowly converging integral
/*! \details This is likely an indication that roundoff error is preventing
 * convergence of the integral.
 */
template<typename IntegralQuantity>
bool AdaptiveIntegrationDivergenceDetector<IntegralQuantity>::isIntegralVerySlowlyConverging(
                                     const IntegralQuantity& original_integral,
                                     const IntegralQuantity& original_error,
                                     const IntegralQuantity& refined_integral,
                                     const IntegralQuantity& refined_error )
{
  // Compute the change in the integral
  IntegralQuantity integral_change =
    fabs(original_integral - refined_integral);


  // Check if the integral change limit has been reached
  const bool integral_change_limit_reached =
    integral_change <= d_integral_change_upper_limit*fabs(refined_integral);

  // Check if the error slow change limit has been reached
  const bool error_slow_change_limit_reached =
    refined_error >= d_error_slow_change_lower_limit*original_error;

  // If both limits were reached simultaneously the integral is very
  // slowly converging
  return (integral_change_limit_reached && error_slow_change_limit_reached);
}

// Check for a strictly divering integral
/*! \details A growing error will only be treated as an indication of a 
 * diverging integral if the iteration limit has been reached.
 */
template<typename IntegralQuantity>
bool
AdaptiveIntegrationDivergenceDetector<IntegralQuantity>::isIntegralDiverging(
                                        const IntegralQuantity& original_error,
                                        const IntegralQuantity& refined_error,
                                        const unsigned number_of_iterations )
{
  // Iteration limit reached - check for a growing error
  if( number_of_iterations >= d_iteration_limit )
    return refined_error > original_error;

  // Iteration limit not reached - ignore error behavior
  else
    return false;          
}
  
} // end Utility namespace

#endif // end UTILITY_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_AdaptiveIntegrationDivergenceDetector_def.hpp
//---------------------------------------------------------------------------//
