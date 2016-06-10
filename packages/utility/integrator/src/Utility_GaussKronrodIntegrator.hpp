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
#include "Utility_ExtrapolatedQuadratureBin.hpp"
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

  //! Constructor
  UnitAwareGaussKronrodIntegrator( const FloatType relative_error_tol,
                                   const IntegralQuantity absolute_error_tol =
                                   IntegralQuantity::zero(),
                                   const size_t subinterval_limit = 1000 );

  //! Destructor
  ~UnitAwareGaussKronrodIntegrator()
  { /* ... */ }

  //! Integrate the function adaptively with BinQueue
  template<int Points, typename Functor>
  void integrateAdaptively( Functor& integrand,
			    ArgQuantity lower_limit,
			    ArgQuantity upper_limit,
			    IntegralQuantity& integral,
			    IntegralQuantity& absolute_error ) const;

  //! Integrate the function with point rule
  template<int Points, typename Functor>
  void integrateWithPointRule( Functor& integrand,
                               ArgQuantity lower_limit,
                               ArgQuantity upper_limit,
                               IntegralQuantity& integral,
                               IntegralQuantity& absolute_error,
                               IntegralQuantity& integral_abs, 
                               IntegralQuantity& integral_asc ) const;

  //! Integrate a function with known integrable singularities adaptively
  template<typename Functor, typename ArrayType>
  void integrateAdaptivelyWynnEpsilon( Functor& integrand,
                                       const ArrayType& points_of_interest,
                                       IntegralQuantity& integral,
                                       IntegralQuantity& absolute_error) const;

protected:

  // The quadrature bin type
  typedef QuadratureBin<ArgQuantity,IntegralQuantity> QuadratureBinType;

  // The extrapolated quadrature bin type
  typedef ExtrapolatedQuadratureBin<ArgQuantity,IntegralQuantity> ExtrapolatedQuadratureBinType;

  // Conduct the first iteration of the adapative integration
  template<int Points, typename Functor>
  bool integrateAdaptivelyInitialIteration( Functor& integrand,
                                            QuadratureBinType& bin ) const;

  // Conduct the other iterations of the adaptive integration
  template<int Points, typename Functor>
  void integrateAdaptivelyIterate( Functor& integrand,
                                   QuadratureBinType& initial_bin,
                                   IntegralQuantity& result,
                                   IntegralQuantity& absolute_error ) const;

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
                             IntegralQuantity integral_abs, 
                             IntegralQuantity integral_asc ) const;

  // Initialize the bins for the adaptive Wynn Epsilon integration
  template<int Points, typename Functor, typename ArrayType>
  bool initializeBinsWynnEpsilon( Functor& integrand,
                                  const ArrayType& points_of_interest,
                                  std::vector<unsigned>& bin_order,
                                  IntegralQuantity& result,
                                  IntegralQuantity& absolute_error ) const;
  
  // Bisect and integrate the given bin interval
  template<int Points, typename Functor, typename Bin>
  void bisectAndIntegrateBinInterval( 
    Functor& integrand, 
    const Bin& bin,
    Bin& bin_1,
    Bin& bin_2,
    IntegralQuantity& bin_1_asc,
    IntegralQuantity& bin_2_asc ) const;

  // Test if subinterval is too small
  template<int Points>
  bool subintervalTooSmall( ArgQuantity lower_limit_1, 
                            ArgQuantity lower_limit_2, 
                            ArgQuantity upper_limit_2 ) const;

  // Check the roundoff error 
  void checkRoundoffError( 
                       const QuadratureBinType& bin, 
                       const QuadratureBinType& bin_1, 
                       const QuadratureBinType& bin_2,    
                       const IntegralQuantity bin_1_asc,
                       const IntegralQuantity bin_2_asc,
                       int& round_off_1,
                       int& round_off_2,
                       const int number_of_iterations ) const;

  // Check the roundoff error 
  void checkRoundoffError( 
                       const ExtrapolatedQuadratureBin<T>& bin, 
                       const ExtrapolatedQuadratureBin<T>& bin_1, 
                       const ExtrapolatedQuadratureBin<T>& bin_2,    
                       const IntegralQuantity bin_1_asc,
                       const IntegralQuanttiy bin_2_asc,
                       int& round_off_1,
                       int& round_off_2,
                       int& round_off_3,
                       const bool extrapolate, 
                       const int number_of_iterations ) const;
 
  // Sort the bin order from highest to lowest error 
  void sortBins( 
        std::vector<int>& bin_order,
        BinArray& bin_array, 
        const ExtrapolatedQuadratureBinType& bin_1,
        const ExtrapolatedQuadratureBinType& bin_2,
        const int number_of_intervals,
        int nr_max ) const;

  // Get the Wynn Epsilon-Algoirithm extrapolated value
  template<typename Array>
  void getWynnEpsilonAlgorithmExtrapolation( 
                                    Array& bin_extrapolated_result, 
                                    Array& last_three_results, 
                                    IntegralType& extrapolated_result, 
                                    IntegralType& extrapolated_error,  
                                    int& number_of_extrapolated_intervals,
                                    int& number_of_extrapolated_calls  ) const;

private:

  // The quadrature bin priority queue
  typedef std::priority_queue<QuadratureBinType> BinQueue;

  // The extrapolated quadrature bin array
  typedef std::vector<ExtrapolatedQuadratureBinType> BinArray;
  
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
