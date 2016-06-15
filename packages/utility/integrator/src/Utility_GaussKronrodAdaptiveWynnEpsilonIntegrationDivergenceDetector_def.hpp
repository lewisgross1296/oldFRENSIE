//---------------------------------------------------------------------------//
//!
//! \file   Utility_GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector_def.hpp
//! \author Alex Robinson
//! \brief  The Gauss-Kronrod adaptive Wynn-Epsilon int. div. detector template
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_GAUSS_KRONROD_ADAPTIVE_WYNN_EPSILON_INTEGRATION_DIVERGENCE_DETECTOR_DEF_HPP
#define UTILITY_GAUSS_KRONROD_ADAPTIVE_WYNN_EPSILON_INTEGRATION_DIVERGENCE_DETECTOR_DEF_HPP

// FRENSIE Includes
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_ContractException.hpp"

namespace Utility{

// Constructor
template<typename IntegralQuantity>
GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector<IntegralQuantity>::GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector(
            const unsigned slowly_converging_extrapolated_integral_count_limit,
            const unsigned slowly_converging_integral_count_limit ,
            const unsigned diverging_integral_count_limit,
            const typename IQT::RawType integral_change_upper_limit,
            const typename IQT::RawType error_slow_change_lower_limit,
            const unsigned iteration_limit )
  : GaussKronrodAdaptiveIntegrationDivergenceDetector(
                                        slowly_converging_integral_count_limit,
                                        diverging_integral_count_limit,
                                        integral_change_upper_limit,
                                        error_slow_change_lower_limit,
                                        interation_limit ),
    d_slowly_converging_extrapolated_integral_counter( 0u ),
    d_slowly_converging_extrapolated_integral_count_limit(
                         slowly_converging_extrapolated_integral_count_limit ),
    d_extrapolated_integrals( false )
    
{
  // Make sure the slowly converging extrapolated integral count limit is valid
  testPrecondition( slowly_converging_extrapolated_integral_count_limit > 0u );
}

// Indicate that integrals are being extrapolated
template<typename IntegralQuantity>
void GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector<IntegralQuantity>::setExtrapolatedIntegrals( const bool extrapolated )
{
  d_extrapolated_integrals = extrapolated;
}

// Check for a diverging integral estimate
template<typename IntegralQuantity>
void GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector<IntegralQuantity>::checkForDivergence(
                                     const IntegralQuantity& original_integral,
                                     const IntegralQuantity& original_error,
                                     const IntegralQuantity& refined_integral,
                                     const IntegralQuantity& refined_error,
                                     const unsigned number_of_iterations )
{
  // Do the basic check
  GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::checkForDivergence(
                                                        original_integral,
                                                        original_error,
                                                        refined_integral,
                                                        refined_error,
                                                        number_of_iterations );

  // Check if the slowly converging extrapolated integral count limit has
  // been reached
  if( d_slowly_converging_extrapolated_integral_counter >=
      d_slowly_converging_extrapolated_integral_count_limit )
  {
    THROW_EXCEPTION( Utility::IntegratorException,
                     "Error: Extremely bad integrand behavior occurs at some "
                     "points of the integration interval! Convergence cannot "
                     "be achieved (" << number of iterations <<
                     " iterations)!" );
  }
}

// Get the slowly converging extrapolated integral count
template<typename IntegralQuantity>
unsigned GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector<IntegralQuantity>::getSlowlyConvergingExtrapolatedIntegralCount() const
{
  return d_slowly_converging_extrapolated_integral_counter;
}

// Get the slowly converging extrapolated integral count limit
template<typename IntegralQuantity>
unsigned GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector<IntegralQuantity>::getSlowlyConvergingExtrapolatedIntegralCountLimit() const
{
  return d_slowly_converging_extrapolated_integral_count_limit;
}

// Increment the slowly converging integral counter
/*! \details If the extrapolation flag has been set, increment the
 * slowly converging extrapolated integral counter. Otherwise, increment
 * the slowly converging integral counter.
 */
template<typename IntegralQuantity>
void GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector<IntegralQuantity>::incrementSlowlyConvergingIntegralCounter()
{
  if( d_extrapolated_integrals )
    ++d_slowly_converging_extrapolated_integral_counter;
  else
  {
    GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::incrementSlowlyConvergingIntegralCounter();
  }
}

// Test if the slowly converging integral count limit has been reached
/*! \details The sum of the extrapolated and normal slowly converging
 * integral counters will be used to conduct the test.
 */
template<typename IntegralQuantity>
bool GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector<IntegralQuantity>::hasSlowlyConvergingIntegralCountLimitBeenReached() const
{
  return this->getSlowlyConvergingIntegralCount() +
    d_slowly_converging_extrapolated_integral_counter >=
    this->getSlowlyConvergingIntegralCountLimit();
}
  
} // end Utility namespace

#endif // end UTILITY_GAUSS_KRONROD_ADAPTIVE_WYNN_EPSILON_INTEGRATION_DIVERGENCE_DETECTOR_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector_def.hpp
//---------------------------------------------------------------------------//
