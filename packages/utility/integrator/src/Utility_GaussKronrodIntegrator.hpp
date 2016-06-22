//---------------------------------------------------------------------------//
//!
//! \file   Utility_GaussKronrodIntegrator.hpp
//! \author Luke Kersting
//! \brief  Gauss-Kronrod integrator
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_GAUSS_KRONROD_INTEGRATOR_HPP
#define UTILITY_GAUSS_KRONROD_INTEGRATOR_HPP

// Std Includes
#include <queue>
#include <vector>

// FRENSIE Includes
#include "Utility_QuadratureBin.hpp"
#include "Utility_WynnEpsilonExtrapolationTable.hpp"
#include "Utility_WynnEpsilonPriorityQueue.hpp"
#include "Utility_UnitTraits.hpp"
#include "Utility_QuantityTraits.hpp"

namespace Utility{

//! The Gauss-Kronrod integrator
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
class UnitAwareGaussKronrodIntegrator
{

private:

  // The argument unit traits typedef
  typedef UnitTraits<ArgUnit> ArgUnitTraits;

  // The integrand unit traits typedef
  typedef UnitTraits<IntegrandUnit> IntegrandUnitTraits;

  // The integral unit traits typedef
  typedef UnitTraits<typename IntegrandUnitTraits::template GetMultipliedUnitType<ArgUnit>::type> IntegralUnitTraits;

  // The argument quantity traits typedef
  typedef QuantityTraits<typename ArgUnitTraits::template GetQuantityType<FloatType>::type> ArgQT;

  // The integrand quantity traits typedef
  typedef QuantityTraits<typename IntegrandUnitTraits::template GetQuantityType<FloatType>::typ>e IntegrandQT;

  // The integral quantity traits typedef
  typedef QuantityTraits<typename IntegralUnitTraits::template GetQuantityType<FloatType>::type> IntegralQT;

  // The unitless quantity traits typedef
  typedef QuantityTraits<FloatType> QT;

public:

  //! The argument quantity type
  typedef typename ArgQT::QuantityType ArgQuantity;

  //! The integrand quantity type
  typedef typename IntegrandQT::QuantityType IntegrandQuantity;

  //! The integral quantity type
  typedef typename IntegralQT::QuantityType IntegralQuantity;

  //! The quadrature bin type
  typedef QuadratureBin<ArgQuantity,IntegralQuantity> QuadratureBinType;

  //! Constructor
  UnitAwareGaussKronrodIntegrator( const FloatType relative_error_tol,
                                   const IntegralQuantity absolute_error_tol =
                                   IntegralQuantity::zero(),
                                   const size_t subinterval_limit = 1000 );

  //! Destructor
  ~UnitAwareGaussKronrodIntegrator()
  { /* ... */ }

  //! Integrate the function with point rule
  template<int Points, typename Functor>
  void integrateWithPointRule( Functor& integrand,
                               QuadratureBinType& bin ) const;

  //! Integrate the function with point rule
  template<int Points, typename Functor>
  void integrateWithPointRule( Functor& integrand,
                               const ArgQuantity lower_limit,
                               const ArgQuantity upper_limit,
                               IntegralQuantity& integral,
                               IntegralQuantity& absolute_error,
                               IntegralQuantity& integral_abs, 
                               IntegralQuantity& integral_asc ) const;

  //! Integrate the function adaptively with point rule
  template<int Points, typename Functor>
  void integrateAdaptively( Functor& integrand,
                            QuadratureBinType& bin ) const;
  
  //! Integrate the function adaptively with point rule
  template<int Points, typename Functor>
  void integrateAdaptively( Functor& integrand,
			    const ArgQuantity lower_limit,
			    const ArgQuantity upper_limit,
			    IntegralQuantity& integral,
			    IntegralQuantity& absolute_error ) const;

  //! Integrate a function with known integrable singularities adaptively
  template<typename Functor, typename ArrayType>
  void integrateAdaptivelyWynnEpsilon( Functor& integrand,
                                       const ArrayType& points_of_interest,
                                       IntegralQuantity& integral,
                                       IntegralQuantity& absolute_error) const;

protected:

  // The quadrature bin priority queue
  typedef std::priority_queue<QuadratureBinType,std::vector<QuadratureBinType>,QuadratureBinType::compareBinErrors> BinQueue;

  // The Wynn-Epsilon quadrature bin priority queue
  typedef WynnEpsilonPriorityQueue<ArgQuantity,IntegralQuantity> WynnEpsilonBinQueue;

  // The Wynn-Epslion extrapolation table type
  typedef WynnEpsilonExtrapolationTable<IntegralQuantity> WynnEpsilonExtrapolationTableType;

  // Evaluate the integrand at the kronrod abscissae
  template<typename Functor>
  IntegralQuantity evaluateIntegrandAtKronrodAbscissae(
                     Functor& integrand,
                     std::vector<IntegrandQuantity>& integrand_values_neg_absc,
                     std::vector<IntegrandQuantity>& integrand_values_pos_absc,
                     IntegrandQuantity& integrand_value_midpoint_absc,
                     const ArgQuantity scale_factor,
                     const ArgQuantity shift_factor,
                     const std::vector<FloatType>& kronrod_abscissae ) const;

  // Calculate the unscaled kronrod integral
  IntegralQuantity calculateUnscaledKronrodIntegral(
               const std::vector<IntegrandQuantity>& integrand_values_neg_absc,
               const std::vector<IntegrandQuantity>& integrand_values_pos_absc,
               const IntegrandQuantity& integrand_value_midpoint_absc,
               const std::vector<FloatType>& kronrod_weights ) const;

  // Calculate the unscaled kronrod integral abs (used for absolute error calc)
  IntegralQuantity calculateUnscaledKronrodIntegralAbs(
               const std::vector<IntegrandQuantity>& integrand_values_neg_absc,
               const std::vector<IntegrandQuantity>& integrand_values_pos_absc,
               const IntegrandQuantity& integrand_value_midpoint_absc,
               const std::vector<FloatType>& kronrod_weights ) const;

  // Calculate the unscaled kronrod integral asc (used for absolute error calc)
  IntegralQuantity calculateUnscaledKronrodIntegralAsc(
               const std::vector<IntegrandQuantity>& integrand_values_neg_absc,
               const std::vector<IntegrandQuantity>& integrand_values_pos_absc,
               const IntegrandQuantity& integrand_value_midpoint_absc,
               const std::vector<FloatType>& kronrod_weights,
               const IntegrandQuantity& unscaled_kronrod_integral ) const;

  // Calculate the unscaled gauss integral 
  IntegralQuantity calculateUnscaledGaussIntegral(
               const std::vector<IntegrandQuantity>& integrand_values_neg_absc,
               const std::vector<IntegrandQuantity>& integrand_values_pos_absc,
               const IntegrandQuantity& integrand_value_midpoint_absc,
               const std::vector<FloatType>& gauss_weights ) const;

  // Rescale absolute error from integration
  void rescaleAbsoluteError( IntegralQuantity& absolute_error, 
                             const IntegralQuantity& integral_abs, 
                             const IntegralQuantity& integral_asc ) const;

  // Conduct the first iteration of the adapative integration
  template<int Points, typename Functor>
  bool integrateAdaptivelyInitialIteration( Functor& integrand,
                                            QuadratureBinType& bin ) const;

  // Conduct the other iterations of the adaptive integration
  template<int Points, typename Functor>
  void integrateAdaptivelyIterate( Functor& integrand,
                                   QuadratureBinType& initial_bin ) const;

  // Test if subinterval is too small
  template<int Points>
  bool subintervalTooSmall( const QuadratureBinType& left_bin,
                            const QuadratureBinType& right_bin ) const;

  // Initialize the bins for the adaptive Wynn Epsilon integration
  template<typename Functor, typename ArrayType>
  bool initializeBinsWynnEpsilon( Functor& integrand,
                                  const ArrayType& points_of_interest,
                                  WynnEpsilonBinQueue& bin_queue,
                                  IntegralQuantity& integral,
                                  IntegralQuantity& absolute_error,
                                  IntegralQuantity& total_absolute_error,
                                  IntegralQuantity& tolerance,
                                  bool& positive_integrand ) const;

  // Conduct the other iterations of the adaptive Wynn Epsilon integration
  template<typename Functor>
  void integrateAdaptivelyWynnEpsilonIterate(
                                     Functor& integrand,
                                     WynnEpsilonBinQueue& bin_queue,
                                     IntegralQuantity& integral,
                                     IntegralQuantity& absolute_error,
                                     IntegralQuantity& total_absolute_error,
                                     const IntegralQuantity& initial_tolerance,
                                     const bool positive_integrand ) const;

private:
  
  // The relative error tolerance
  FloatType d_relative_error_tol;

  // The absolute error tolerance
  IntegralType d_absolute_error_tol;

  // The subinterval limit
  size_t d_subinterval_limit;
};

//! The GaussKronrodIntegrator (unit-agnostic)
template<typename FloatType> using GaussKronrodIntegrator =
  UnitAwareGaussKronrodIntegrator<FloatType,void,void>;

} // end Utility namespace

//---------------------------------------------------------------------------//
// Template Includes
//---------------------------------------------------------------------------//

#include "Utility_GaussKronrodIntegrator_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_GUASS_KRONROD_INTEGRATOR_HPP

//---------------------------------------------------------------------------//
// end Utility_GaussKronrodIntegrator.hpp
//---------------------------------------------------------------------------//
