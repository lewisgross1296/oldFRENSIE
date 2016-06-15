//---------------------------------------------------------------------------//
//!
//! \file   Utility_GaussKronrodAdaptiveIntegrationDivergenceDetector.hpp
//! \author Alex Robinson
//! \brief  The Gauss-Kronrod adaptive integration divergence detector decl.
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_GAUSS_KRONROD_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_HPP
#define UTILITY_GAUSS_KRONROD_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_HPP

// FRENSIE Includes
#include "Utility_QuantityTraits.hpp"
#include "Utility_AdaptiveIntegrationDivergenceDetector.hpp"

namespace Utility{
  
//! The Gauss-Kronrod adaptive integration detector class
template<typename IntegralQuantity>
class GaussKronrodAdaptiveIntegrationDivergenceDetector : public AdaptiveIntegrationDivergenceDetector<IntegralQuantity>
{

protected:

  // Quantity traits for the IntegralQuantity
  typedef AdaptiveIntegrationDivergenceDetector<IntegralQuantity>::IQT IQT;

  // Quantity traits for the IntegralQuantity RawType
  typedef AdaptiveIntegrationDivergenceDetector<IntegralQuantity>::QT QT;

public:

  //! Constructor
  GaussKronrodAdaptiveIntegrationDivergenceDetector(
                   const unsigned slowly_converging_integral_count_limit = 6u,
                   const unsigned diverging_integral_count_limit = 20u,
                   const typename IQT::RawType integral_change_upper_limit =
                   QT::one()/100000,
                   const typename IQT::RawType error_slow_change_lower_limit =
                   QT::one()*99/100;
                   const unsigned iteration_limit = 10u );

  //! Destructor
  virtual ~GaussKronrodAdaptiveIntegrationDivergenceDetector()
  { /* ... */ }

  //! Check for a diverging integral estimate
  virtual void checkForDivergence( const IntegralQuantity& original_integral,
                                   const IntegralQuantity& original_error,
                                   const IntegralQuantity& refined_integral,
                                   const IntegralQuantity& refined_error,
                                   const unsigned number_of_iterations );

  //! Get the slowly converging integral count
  unsigned getSlowlyConvergingIntegralCount() const;

  //! Get the diverging integral count
  unsigned getDivergingIntegralCount() const;

  //! Get the slowly converging integral count limit
  unsigned getSlowlyConvergingIntegralCountLimit() const;

  //! Get the diverging integral count limit
  unsigned getDivergingIntegralCountLimit() const;

protected:

  //! Increment the slowly converging integral counter
  virtual void incrementSlowlyConvergingIntegralCounter();

  //! Increment the diverging integral counter
  virtual void incrementDivergingIntegralCounter();

  //! Test if the slowly converging integral count limit has been reached
  virtual bool hasSlowlyConvergingIntegralCountLimitBeenReached() const;

  //! Test if the diverging integral count limit has been reached
  virtual bool hasDivergingIntegralCountLimitBeenReached() const;

private:

  // The slowly converging integral counter
  unsigned d_slowly_converging_integral_counter;

  // The diverging integral counter
  unsigned d_diverging_integral_counter;

  // The slowly converging integral count limit
  unsigned d_slowly_converging_integral_count_limit;

  // The diverging integral count limit
  unsigned d_diverging_integral_count_limit;
};

} // end Utility namespace

//---------------------------------------------------------------------------//
// Template Includes
//---------------------------------------------------------------------------//

#include "Utility_GaussKronrodAdaptiveIntegrationDivergenceDetector_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_GAUSS_KRONROD_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_HPP

//---------------------------------------------------------------------------//
// end Utility_GaussKronrodAdaptiveIntegrationDivergenceDetector.hpp
//---------------------------------------------------------------------------//
