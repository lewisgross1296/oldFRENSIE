//---------------------------------------------------------------------------//
//!
//! \file   Utility_AdaptiveIntegrationDivergenceDetector.hpp
//! \author Alex Robinson
//! \brief  The adaptive integration divergence detector base class declaration
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_HPP
#define UTILITY_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_HPP

// FRENSIE Includes
#include "Utility_QuadartureBin.hpp"
#include "Utility_QuantityTraits.hpp"

namespace Utility{

//! The adaptive integration divergence detector base class
template<typename IntegralQuantity>
class AdaptiveIntegrationDivergenceDetector
{

protected:

  // Quantity traits for the IntegralQuantity
  typedef QuantityTrait<IntegralQuantity> IQT;

  // Quantity traits for the IntegralQuantity RawType
  typedef QuantityTrait<typename IQT::RawType> QT;

public:

  //! Constructor
  AdaptiveIntegrationDivergenceDetector(
                     const typename IQT::RawType integral_change_upper_limit
                     const typename IQT::RawType error_slow_change_lower_limit,
                     const unsigned iteration_limit );

  //! Destructor
  virtual ~AdaptiveIntegrationDivergenceDetector()
  { /* ... */ }

  //! Check for a diverging integral estimate
  virtual void checkForDivergence( const IntegralQuantity& original_integral,
                                   const IntegralQuantity& original_error,
                                   const IntegralQuantity& refined_integral,
                                   const IntegralQuantity& refined_error,
                                   const unsigned number_of_iterations ) = 0;

  //! Check for a diverging integral estimator
  template<typename ArgQuantity>
  void checkForDivergence(
             const QuadratureBin<ArgQuantity,IntegralQuantity>& original_bin,
             const QuadratureBin<ArgQuantity,IntegralQuantity>& left_half_bin,
             const QuadratureBin<ArgQuantity,IntegralQuantity>& right_half_bin,
             const unsigned number_of_iterations );
                           

protected:

  //! Check for a very slowly converging integral 
  bool isIntegralVerySlowlyConverging(
                                     const IntegralQuantity& original_integral,
                                     const IntegralQuantity& original_error,
                                     const IntegralQuantity& refined_integral,
                                     const IntegralQuantity& refined_error );

  //! Check for a strictly divering integral
  bool isIntegralDiverging( const IntegralQuantity& original_error,
                            const IntegralQuantity& refined_error,
                            const unsigned number_of_iterations );
}

private:

  // The integral change upper limit
  // (used to test for a very slowly converging integral)
  typename IQT::RawType d_integral_change_upper_limit;

  // The error slow change lower limit
  // (used to test for a very slowly converging integral)
  typename IQT::RawType d_error_slow_change_lower_limit;

  // The iteration limit
  // (when this limit is reached a growing absolute error will be treated as
  //  an indication of a diverging integral)
  unsigned d_iteration_limit;
};

} // end Utility namespace

//---------------------------------------------------------------------------//
// Template Includes.
//---------------------------------------------------------------------------//

#include "Utility_AdaptiveIntegrationDivergenceDetector_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_ADAPTIVE_INTEGRATION_DIVERGENCE_DETECTOR_HPP

//---------------------------------------------------------------------------//
// end Utility_AdaptiveIntegrationDivergenceDetector.hpp
//---------------------------------------------------------------------------//
