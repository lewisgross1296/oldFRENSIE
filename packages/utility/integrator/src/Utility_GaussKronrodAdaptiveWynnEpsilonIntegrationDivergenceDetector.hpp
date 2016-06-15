//---------------------------------------------------------------------------//
//!
//! \file   Utility_GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector.hpp
//! \author Alex Robinson
//! \brief  The Gauss-Kronrod adaptive Wynn-Epsilon int. div. detector decl.
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_GAUSS_KRONROD_ADAPTIVE_WYNN_EPSILON_INTEGRATION_DIVERGENCE_DETECTOR_HPP
#define UTILITY_GAUSS_KRONROD_ADAPTIVE_WYNN_EPSILON_INTEGRATION_DIVERGENCE_DETECTOR_HPP

// FRENSIE Includes
#include "Utility_QuantityTraits.hpp"
#include "Utility_GaussKronrodAdaptiveIntegrationDivergenceDetector.hpp"

namespace Utility{

//! The Gauss-Kronrod adaptive Wynn-Epsilon integration detector class
template<typename IntegralQuantity>
class GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector : public GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>
{

protected:

  // Quantity traits for the IntegralQuantity
  typedef GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::IQT IQT;

  // Quantity traits for the IntegralQuantity RawType
  typedef GaussKronrodAdaptiveIntegrationDivergenceDetector<IntegralQuantity>::QT QT;

public:

  //! Constructor
  GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector(
       const unsigned slowly_converging_extrapolated_integral_count_limit = 5u,
       const unsigned slowly_converging_integral_count_limit = 10u,
       const unsigned diverging_integral_count_limit = 20u,
       const typename IQT::RawType integral_change_upper_limit =
       QT::one()/100000,
       const typename IQT::RawType error_slow_change_lower_limit =
       QT::one()*99/100;
       const unsigned iteration_limit = 10u );

  //! Destructor
  ~GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector()
  { /* ... */ }

  //! Indicate that integrals are being extrapolated
  void setExtrapolatedIntegrals( const bool extrapolated );

  //! Check for a diverging integral estimate
  void checkForDivergence( const IntegralQuantity& original_integral,
                           const IntegralQuantity& original_error,
                           const IntegralQuantity& refined_integral,
                           const IntegralQuantity& refined_error,
                           const unsigned number_of_iterations );

  //! Get the slowly converging extrapolated integral count
  unsigned getSlowlyConvergingExtrapolatedIntegralCount() const;

  //! Get the slowly converging extrapolated integral count limit
  unsigned getSlowlyConvergingExtrapolatedIntegralCountLimit() const;

protected:

  //! Increment the slowly converging integral counter
  void incrementSlowlyConvergingIntegralCounter();

  //! Test if the slowly converging integral count limit has been reached
  bool hasSlowlyConvergingIntegralCountLimitBeenReached() const;

private:

  // The slowly converging extrapolated integral counter
  unsigned d_slowly_converging_extrapolated_integral_counter;

  // The slowly converging extraplated integral count limit
  unsigned d_slowly_converging_extrapolated_integral_count_limit;

  // Treat the integrals as extrapolated (if True)
  bool d_extrapolated_integrals;

} // end Utility namespace

#endif // end UTILITY_GAUSS_KRONROD_ADAPTIVE_WYNN_EPSILON_INTEGRATION_DIVERGENCE_DETECTOR_HPP

//---------------------------------------------------------------------------//
// end Utility_GaussKronrodAdaptiveWynnEpsilonIntegrationDivergenceDetector.hpp
//---------------------------------------------------------------------------//
