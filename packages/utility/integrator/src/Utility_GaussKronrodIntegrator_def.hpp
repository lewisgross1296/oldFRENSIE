//---------------------------------------------------------------------------//
//!
//! \file   Utility_GaussKronrodIntegrator_def.hpp
//! \author Luke Kersting
//! \brief  Gauss-Kronrod integrator
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_GAUSS_KRONROD_INTEGRATOR_DEF_HPP
#define UTILITY_GAUSS_KRONROD_INTEGRATOR_DEF_HPP

// Std Includes
#include <limits>
#include <algorithm>

// FRENSIE Includes
#include "Utility_ContractException.hpp"
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_ExceptionCatchMacros.hpp"
#include "Utility_IntegratorException.hpp"
#include "Utility_SortAlgorithms.hpp"
#include "Utility_QuantityTraits.hpp"
#include "Utility_GaussKronrodQuadratureSetTraits.hpp"

namespace Utility{

// Constructor
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::UnitAwareGaussKronrodIntegrator( 
    const FloatType relative_error_tol,
    const IntegralQuantity absolute_error_tol,
    const size_t subinterval_limit )
  : d_relative_error_tol( relative_error_tol ),
    d_absolute_error_tol( absolute_error_tol ),
    d_subinterval_limit( subinterval_limit )
{
  // Make sure the floating point type is valid
  testStaticPrecondition( (boost::is_floating_point<FloatType>::value) );
  // Make sure the error tolerances are valid
  testPrecondition( relative_error_tol >= QT::zero() );
  testPrecondition( absolute_error_tol >= IntegralQT::zero() );
  // Make sure the subinterval limit is valid
  testPrecondition( subinterval_limit > 0 );

  // Check if the tolerance can be achieved
  if( d_absolute_error_tol <= QT::zero() )
  {
    THROW_EXCEPTION( Utility::IntegratorException,
                     "Error: Convergence will not be possible with "
                     "the given absolute error tolerance!" );
  }

  if( d_relative_error_tol < 50*QT::epsilon() )
  {
    THROW_EXCEPTION( Utility::IntegratorException,
                     "Error: Convergence will not be possible with "
                     "the given relative error tolerance!" );
  }
}

// Integrate the function with given Gauss-Kronrod point rule
/*! \details See the overloaded method for more details.
 */
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<int Points, typename Functor>
inline void UnitAwareGaussKronrodIntegrator<T>::integrateWithPointRule(
                                         Functor& integrand, 
                                         ArgQuantity lower_limit, 
                                         ArgQuantity upper_limit,
                                         IntegralQuantity& integral,
                                         IntegralQuantity& absolute_error,
                                         IntegralQuantity& integral_abs, 
                                         IntegralQuantity& integral_asc ) const
{
  // Make sure the integration limits are valid
  testPrecondition( lower_limit <= upper_limit );

  // Create the quadrature bin for the limits of interest
  QuadratureBinType bin( lower_limit, upper_limit );

  this->integrateWithPointRule<Points>( integrand, bin );

  integral = bin.getIntegral();
  absolute_error = bin.getAbsoluteError();
  integral_abs = bin.getIntegralAbs();
  integral_asc = bin.getIntegralAsc();
}

// Integrate the function with given Gauss-Kronrod point rule
/*! \details Functor must have operator()( double ) defined. This function
 * applies the specified integration rule (Points) to estimate 
 * the integral of the integrand over [lower_limit,upper_limit]. 
 * Valid Gauss-Kronrod rules are 15, 21, 31, 41, 51 and 61. 
 * Higher-order rules give better accuracy for smooth functions, 
 * while lower-order rules save time when the function contains local 
 * difficulties, such as discontinuities. On each iteration
 * the adaptive integration strategy bisects the interval with the largest
 * error estimate. See the qag function details in the quadpack documentation.
 */
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<int Points, typename Functor>
void UnitAwareGaussKronrodIntegrator<T>::integrateWithPointRule(
                                                 Functor& integrand, 
                                                 QuadratureBinType& bin ) const
{
  // Make sure the integration limits are valid
  testPrecondition( bin.getLowerLimit() <= bin.getUpperLimit() );
  
  if( bin.getLowerLimit() < bin.getUpperLimit() )
  {
    // The Gauss integral
    IntegralQuantity gauss_integral = IntegralQT::zero();

    // Calculate the scale factor
    ArgQuantity scale_factor = (upper_limit - lower_limit)/2;

    // Calculate the shift factor
    ArgQuantity shift_factor = (upper_limit + lower_limit)/2;

    // Calculate the unscaled integrals, then scale
    {
      // Get the quadrature set traits for the requested point rule
      typedef GaussKronrodQuadratureSetTraits<FloatType> GKQST;
      
      // Evaluate the integrand at the kronrod abscissae
      std::vector<IntegrandQuantity> integrand_values_neg_absc;
      std::vector<IntegrandQuantity> integrand_values_pos_absc;
      IntegrandQuantity integrand_value_midpoint_absc;

      this->evaluateIntegrandAtKronrodAbscissae( integrand,
                                                 integrand_values_neg_absc,
                                                 integrand_values_pos_absc,
                                                 integrand_value_midpoint_absc,
                                                 scale_factor,
                                                 shift_factor,
                                                 GKQST::getKronrodAbscissae());

      // Calculate the unscaled kronrod integral
      IntegrandQuantity unscaled_kronrod_integral =
        this->calculateUnscaledKronrodIntegral( integrand_values_neg_absc,
                                                integrand_values_pos_absc,
                                                integrand_value_midpoint_absc,
                                                GKQST::getKronrodWeights() );

      // Calculate the unscaled kronrod integral abs
      IntegrandQuantity unscaled_kronrod_integral_abs =
        this->calculateUnscaledKronrodIntegralAbs(
                                                integrand_values_neg_absc,
                                                integrand_values_pos_absc,
                                                integrand_values_midpoint_absc,
                                                GKQST::getKronrodWeights() );

      // Calculate the unscaled kronrod integral asc
      IntegrandQuantity unscaled_kronrod_integral_asc =
        this->calculateUnscaledKronrodIntegralAsc(
                                                integrand_values_neg_absc,
                                                integrand_values_pos_absc,
                                                integrand_values_midpoint_absc,
                                                GKQST::getKronrodWeights(),
                                                unscaled_kronrod_integral );   

      // Calculate the unscaled gauss integral
      IntegrandQuantity unscaled_gauss_integral =
        this->calculateUnscaledGaussIntegral( integrand_values_neg_absc,
                                              integrand_values_pos_absc,
                                              integrand_value_midpoint_absc,
                                              GKQST::getGaussWeights() );

      // Scale the integrals
      bin.setIntegral( unscaled_kronrod_integral*scale_factor );
      bin.setIntegralAbs( unscaled_kronrod_integral_abs*fabs( scale_factor ) );
      bin.setIntegralAsc( unscaled_kronrod_integral_asc*fabs( scale_factor ) );
      
      gauss_integral = unscaled_gauss_integral*scale_factor;
    }

    // Estimate error in integral 
    bin.setAbsoluteError( fabs( bin.getIntegral() - gauss_integral ) );
    
    this->rescaleAbsoluteError( bin.getAbsoluteError(),
                                bin.getIntegralAbs(),
                                bin.getIntegralAsc() );
  }
  else if( lower_limit == upper_limit )
  {
    bin.setIntegral( IntegralQT::zero() );
    bin.setIntegralAbs( IntegralQT::zero() );
    bin.setIntegralAsc( IntegralQT::zero() );
    bin.setAbsoluteError( IntegralQT::zero() );
  }
  else // Invalid limits
  {
    THROW_EXCEPTION( Utility::IntegratorException,
		     "Error: Invalid integration limits ( "
                     << lower_limit << " !< " << upper_limit << ")." );
  }
}

// Integrate the function adaptively
/*! \details See overloaded method for more details.
 */
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<int Points, typename Functor>
inline void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::integrateAdaptively(
                                       Functor& integrand, 
                                       ArgQuantity lower_limit, 
                                       ArgQuantity upper_limit,
                                       IntegralQuantity& integral,
				       IntegralQuantity& absolute_error ) const
{
  // Make sure the integration limits are valid
  testPrecondition( lower_limit <= upper_limit );
  
  // Create the quadrature bin for the limits of interest
  QuadratureBinType bin( lower_limit, upper_limit );

  this->integrateAdaptively<Points>( integrand, bin );

  integral = bin.getIntegral();
  absolute_error = bin.getAbsoluteError();
}

// Integrate the function adaptively 
/*! \details Functor must have operator()( double ) defined. This function
 * applies the specified integration rule (Points) adaptively until an 
 * estimate of the integral of the integrand over [lower_limit,upper_limit] is
 * achieved within the desired tolerances. Valid Gauss-Kronrod rules are
 * 15, 21, 31, 41, 51 and 61. Higher-order rules give better accuracy for
 * smooth functions, while lower-order rules save time when the function
 * contains local difficulties, such as discontinuities. On each iteration
 * the adaptive integration strategy bisects the interval with the largest
 * error estimate. See the qag function details in the quadpack documentation.
 */
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<int Points, typename Functor>
void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::integrateAdaptively(
                                                 Functor& integrand, 
                                                 QuadratureBinType& bin ) const
{
  // Make sure the integration limits are valid
  testPrecondition( bin.getLowerLimit() <= bin.getUpperLimit() );
  
  bool converged;

  try{
    converged = this->integrateAdaptivelyInitialIteration( integrand, bin );
  }
  EXCEPTION_CATCH_RETHROW( Utility::IntegratorException,
                           "Error: The initial iteration of the adaptive "
                           "integration failed!" );
  
  // Check for fast convergence
  if( converged )
    return;

  // Check if the subinterval limit is insufficient for convergence
  else if( d_subinterval_limit == 1 )
  {
    THROW_EXCEPTION( Utility::IntegratorException,
                     "Error: a maximum of one iteration was insufficient!" );
  }
  
  // Iterate until convergence (or reaching subinterval limit)
  else
  {
    try{
      this->integrateAdaptivelyIterate( integrand, bin );
    }
    EXCEPTION_CATCH_RETHROW( Utility::IntegratorException,
                             "Error: The integration failed during the "
                             "subinterval iterations!" );
  }
}

// Integrate a function with integrable singularities adaptively
/*! \details Functor must have operator()( double ) defined. This function
 * applies the Gauss-Kronrod 21-point integration rule adaptively until an
 * estimate of the integral of the integrand over [lower_limit,upper_limit]
 * is achieved within the desired tolerances. The results are extrapolated 
 * using the Wynn epsilon-algorithm, which accelerates the convergence of the 
 * integral in the presence of discontinuities and integrable singularities. 
 * See QAGS algorithm details in the GNU Scientific Library documentation.
 */
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<typename Functor, typename ArrayType>
void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::integrateAdaptivelyWynnEpsilon( 
                                       Functor& integrand,
                                       const ArrayType& points_of_interest,
                                       IntegralQuantity& integral,
                                       IntegralQuantity& absolute_error ) const
{
  // Make sure that the ArrayType is valid
  testStaticPrecondition( (boost::is_same<ArrayType::value_type,ArgQuantity>::value) );
  // Check that the points of interest are in ascending order
  testPrecondition( Sort::isSortedAscending( points_of_interest.begin(),
                                             points_of_interest.end(),
                                             true ) );
  // Check that the number of points don't exceed the subinterval limit
  testPrecondition( points_of_interest.size() < d_subinterval_limit );
  // Check that there are at least two points
  testPrecondition( points_of_interest.size() >= 2 );

  // The bin array (corresponding to points of interest and subintervals)
  BinArray bin_array;

  // The total absolute error (>= absolute error)
  IntegralQuantity total_absolute_error = IntegralQT::zero();
  
  // Used for convergence testing
  int ksgn;
  
  // Initialize the bins of the integral and check for fast convergence
  const bool converged;

  try{
    converged = this->initializeBinsWynnEpsilon( integrand,
                                                 points_of_interest,
                                                 bin_array,
                                                 integral,
                                                 absolute_error,
                                                 total_absolute_error,
                                                 ksgn );
  }
  EXCEPTION_CATCH_RETHROW( Utility::IntegratorException,
                           "Error: Could not initialize the intervals of "
                           "the adaptive Wynn Epsilon integration!" );
  
  // The integral has already converged - finish
  if( converged )
    return;
  
  // Check for an insufficient number of iterations
  else if( d_subinterval_limit == 1 )
  {
    THROW_EXCEPTION( Utility::IntegratorException,
                     "Error: A maximum of one iteration was insufficient for "
                     "convergence within the tolerance!" );
  }
  
  // Iterate until convergence (or reaching subinterval limit)
  else
  {
    try{
      this->integrateAdaptivelyWynnEpsilonIterate( integrand,
                                                   bin_array,
                                                   integral,
                                                   absolute_error,
                                                   total_absolute_error,
                                                   points_of_interest.size()-1,
                                                   ksgn );
    }
    EXCEPTION_CATCH_RETHROW( Utility::IntegratorException,
                             "Error: The integration failed during the "
                             "subinterval iterations!" );
  } 
}

// Evaluate the integrand at the kronrod abscissae
// Note: Only the positive abscissae (including 0.0) should be passed to
//       this method.
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<typename Functor>
IntegralQuantity UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::evaluateIntegrandAtKronrodAbscissae(
                     Functor& integrand,
                     std::vector<IntegrandQuantity>& integrand_values_neg_absc,
                     std::vector<IntegrandQuantity>& integrand_values_pos_absc,
                     IntegrandQuantity& integrand_value_midpoint_absc,
                     const ArgQuantity scale_factor,
                     const ArgQuantity shift_factor,
                     const std::vector<FloatType>& kronrod_abscissae ) const
{
  // Make sure the kronrod abscissae are valid
  testPrecondition( kronrod_abscissae.size() > 0 );
  testPrecondition( kronrod_abscissae.back() == (FloatType)0.0 );

  // Get number of positive Kronrod abscissae (not including 0.0)
  unsigned number_of_pos_abscissae = kronrod_abscissae.size()-1;
  
  // Resize the integrand values arrays
  integrand_values_neg_absc.resize( number_of_pos_abscissae );
  integrand_values_pos_absc.resize( number_of_pos_abscissae );

  // Evaluate the integrand at the pos/neg abscissae (not 0.0)
  for( int j = 0; j < number_of_pos_abscissae; j++ )
  {
    integrand_values_neg_absc[j] =
      integrand( -scale_factor*kronrod_abscissae[j] + shift_factor );

    integrand_values_pos_absc[j] =
      integrand( scale_factor*kronrod_abscissae[j] + shift_factor );    
  }

  // Evaluate the integrand at the midpoint abscissa
  integrand_value_midpoint_absc =
    integrand( scale_factor*kronrod_abscissae.back() + shift_factor );
}

// Calculate the unscaled kronrod integral
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
IntegralQuantity UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::calculateUnscaledKronrodIntegral(
               const std::vector<IntegrandQuantity>& integrand_values_neg_absc,
               const std::vector<IntegrandQuantity>& integrand_values_pos_absc,
               const IntegrandQuantity& integrand_value_midpoint_absc,
               const std::vector<FloatType>& kronrod_weights ) const
{
  // Make sure the integrand values and kronrod weights are valid
  testPrecondition( integrand_values_neg_absc.size() ==
                    integrand_values_pos_absc.size() );
  testPrecondition( integrand_values_neg_absc.size() + 1 ==
                    kronrod_weights.size() );

  // The unscaled kronrod integral
  IntegrandQuantity unscaled_kronrod_integral = IntegrandQT::zero();\

  
  unsigned number_of_weights = integrand_values_neg_absc.size()+1;
  
  // Estimate integral using all but last (midpoint) abscissa
  for( int j = 0; j < number_of_weights-1; j++ )
  {  
    // Both +/- abscissa_j use w_j
    unscaled_kronrod_integral += kronrod_weights[j]*
      (integrand_values_neg_absc[j] + integrand_values_pos_absc[j]);
  }
    
  // Update integral with midpoint abscissa value
  unscaled_kronrod_integral +=
    integrand_value_midpoint_absc*kronrod_weights.back();

  return unscaled_kronrod_integral;
}

// Calculate the unscaled kronrod integral abs (used for absolute error calc)
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
IntegralQuantity UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::calculateUnscaledKronrodIntegralAbs(
               const std::vector<IntegrandQuantity>& integrand_values_neg_absc,
               const std::vector<IntegrandQuantity>& integrand_values_pos_absc,
               const IntegrandQuantity& integrand_value_midpoint_absc,
               const std::vector<FloatType>& kronrod_weights ) const
{
  // Make sure the integrand values and kronrod weights are valid
  testPrecondition( integrand_values_neg_absc.size() ==
                    integrand_values_pos_absc.size() );
  testPrecondition( integrand_values_neg_absc.size() + 1 ==
                    kronrod_weights.size() );
  
  // The unscaled kronrod integral
  IntegrandQuantity unscaled_kronrod_integral_abs = IntegrandQT::zero();

  unsigned number_of_weights = integrand_values_neg_absc.size()+1;

  // Estimate integral abs using all but last (midpoint) abscissa
  for( int j = 0; j < number_of_weights-1; j++ )
  {  
    // Both +/- abscissa_j use w_j
    unscaled_kronrod_integral_abs += kronrod_weights[j]*
      (fabs( integrand_values_neg_absc[j] ) +
       fabs( integrand_values_pos_absc[j] ));
  }

  // Update integral abs with midpoint abscissa value
  unscaled_kronrod_integral_abs +=
    fabs( integrand_value_midpoint_absc*kronrod_weights.back() );

  return unscaled_kronrod_integral_abs;
}

// Calculate the unscaled kronrod integral asc (used for absolute error calc)
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
IntegralQuantity UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::calculateUnscaledKronrodIntegralAsc(
               const std::vector<IntegrandQuantity>& integrand_values_neg_absc,
               const std::vector<IntegrandQuantity>& integrand_values_pos_absc,
               const IntegrandQuantity& integrand_value_midpoint_absc,
               const std::vector<FloatType>& kronrod_weights,
               const IntegrandQuantity& unscaled_kronrod_integral ) const
{
  // Make sure the integrand values and kronrod weights are valid
  testPrecondition( integrand_values_neg_absc.size() ==
                    integrand_values_pos_absc.size() );
  testPrecondition( integrand_values_neg_absc.size() + 1 ==
                    kronrod_weights.size() );

  // Calculate the mean kronrod integral
  IntegrandQuantity mean_unscaled_kronrod_integral =
    unscaled_kronrod_integral/2;
    
  // Estimate the integral asc for all but the last weight
  IntegrandQuantity unscaled_kronrod_integral_asc = IntegrandQT::zero();

  unsigned number_of_weights = integrand_values_neg_absc.size()+1;
  
  for( int j = 0; j < number_of_weights - 1; j++ )
  {  
    unscaled_kronrod_integral_asc += kronrod_weights[j]*
      (fabs( integrand_values_neg_absc[j] - mean_unscaled_kronrod_integral ) +
       fabs( integrand_values_pos_absc[j] - mean_unscaled_kronrod_integral ));
  };

  // Estimate the integral asc for the last weight
  unscaled_kronrod_integral_asc += kronrod_weights.back()*
    fabs( integrand_value_midpoint_absc - mean_unscaled_kronrod_integral );

  return unscaled_kronrod_integral_asc;
}
 

// Calculate the unscaled gauss integral using the evaluated integrand values
// Note: The Gauss quadrature set use every other Kronrod abscissa - there
// is no need to re-evaluate the Functor
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
IntegralQuantity UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::calculateUnscaledGaussIntegral(
               const std::vector<IntegrandQuantity>& integrand_values_neg_absc,
               const std::vector<IntegrandQuantity>& integrand_values_pos_absc,
               const IntegrandQuantity& integrand_value_midpoint_absc,
               const std::vector<FloatType>& gauss_weights ) const
{
  // Make sure the integrand values and gauss weights are valid
  testPrecondition( integrand_values_neg_absc.size() ==
                    integrand_values_pos_absc.size() );
  testPrecondition( (integrand_values_neg_absc.size()+1)/2 ==
                    gauss_weights.size() );
  
  IntegrandQuantity unscaled_gauss_integral = IntegrandQT::zero();

  unsigned number_of_kronrod_weights = integrand_values_neg_absc.size()+1;

  unsigned number_of_gauss_weights = number_of_kronrod_weights/2;

  // Estimate integral using all but last (midpoint) abscissa
  for( int j = 0; j < number_of_gauss_weights-1; j++ )
  {
    // Every other Kronrod abscissa is also a Gauss abscissa
    int kronrod_j = j*2 + 1;

    unscaled_gauss_integral += gauss_weights[j]*
      (integrand_values_neg_absc[kronrod_j] +
       integrand_values_pos_absc[kronrod_j]);
  };

  // Update Gauss estimate with last weight if needed (even # of abscissae)
  if( number_of_kronrod_weights % 2 == 0 )
  {
    unscaled_gauss_integral += integrand_value_midpoint_absc*
      gauss_weights.back();
  }

  return unscaled_gauss_integral;
}

// Rescale absolute error from integration
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::rescaleAbsoluteError( 
                                   IntegralQuantity& absolute_error, 
                                   const IntegralQuantity& integral_abs, 
                                   const IntegralQuantity& integral_asc ) const
{
  if( integral_asc != IntegralQT::zero() &&
      absolute_error != IntegralQT::zero() )
  {
    FloatType scale = static_cast<FloatType>(200)*absolute_error/integral_asc;

    if( scale < static_cast<FloatType>(1) )
      absolute_error = integral_asc*pow( scale, 3/(FloatType)2 );
    else
      absolute_error = integral_asc; 
  }

  IntegralQuantity cutoff = IntegralQT::min()/ 
    (static_cast<FloatType>(50)*std::numeric_limits<FloatType>::epsilon())
  
  if( integral_abs > cutoff )
  {
    IntegralQuantity min_error = integral_abs*
      static_cast<FloatType>(50)*std::numeric_limits<FloatType>::epsilon();
      
    if( min_error > absolute_error ) 
      absolute_error = min_error;  
  }
}

// Conduct the first iteration of the adapative integration
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<int Points, typename Functor>
bool UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::integrateAdaptivelyInitialIteration(
						 Functor& integrand, 
						 QuadratureBinType& bin ) const
{
  // Perform an integration using the requested point rule over the interval
  // of interest (specified by the bin)
  this->integrateWithPointRule<Points>( integrand, bin );

  // Compute the convergence tolerance
  IntegralQuantity tolerance =
    std::max( d_absolute_error_tol,
              d_relative_error_tol*fabs(bin.getIntegral()) );

  // Need IEEE rounding here to match original quadpack behavior
  IntegralQuantity round_off_limit =
    50*std::numeric_limits<FloatType>::epsilon()*bin.getIntegralAbs();

  // Test if the desired tolerance can be reached
  TEST_FOR_EXCEPTION( bin.getAbsoluteError() <= round_off_limit &&
                      bin.getAbsoluteError() > tolerance, 
                      Utility::IntegratorException,
                      "Error: Cannot reach tolerance because of roundoff "
                      "error during first iteration!" );

  // The integral has converged
  bool converged;
  
  if( (bin.getAbsoluteError() <= tolerance &&
       bin.getAbsoluteError() != bin.getIntegralAsc()) || 
      bin.getAbsoluteError() == IntegralQT::zero() ) 
    converged = true;
  
  // The integral has not converged - more iterations will be required
  else
    converged = false;

  // Reset the integral abs and integral asc (they will no longer be needed)
  bin.setIntegralAbs( IntegralQT::zero() );
  bin.setIntegralAsc( IntegralQT::zero() );
}

// Conduct the other iterations of the adaptive integration
// Note: The integral and absolute error will be stored in the initial
// quadrature bin.
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<int Points, typename Functor>
void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::integrateAdaptivelyIterate(
                                         Functor& integrand,
                                         QuadratureBinType& initial_bin ) const
{
  // Create a copy of the initial bin to add to the priority queue (the
  // initial bin will be used to store the updated integral and absolute
  // error, so it should not be placed in the priority queue).
  BinQueue bin_queue;
  
  {
    QuadratureBinType initial_bin_copy( initial_bin );
  
    // Initialize the bin priority queue with the initial bin
    bin_queue.push( initial_bin_copy );
  }

  // These test counters will be used to determine if the desired tolerance
  // can be acheived
  unsigned round_off_failed_test_counter_a = 0u;
  unsigned round_off_failed_test_counter_b = 0u;

  // Integrate problematic subintervals until convergence is reached or
  // the subinterval limit is reached (the initial subinterval has been done)
  for( int subinterval = 1; subinterval < d_subinterval_limit; subinterval++ )
  {
    // The bin with the highest error (used for the current subinterval)
    QuadratureBinType problem_bin = bin_queue.top();
    bin_queue.pop();

    // The left and right subinterval bin for the current problem bin
    QuadratureBinType left_half_problem_bin, right_half_problem_bin;

    // Bisect the problem bin and calculate the left and right subinterval
    // bin integrals
    try{
      this->bisectAndIntegrateBinInterval<Points>( integrand, 
                                                   problem_bin,
                                                   left_half_problem_bin,
                                                   right_half_problem_bin );
    }
    EXCEPTION_CATCH_RETHROW( Utility::IntegratorException,
                             "Error: Could not bisect and integrate the "
                             "subinterval (" << problem_bin.getLowerLimit() <<
                             "," << problem_bin.getUpperLimit() << ")!" );

    // Check if the subinterval has gotten too small
    if( this->subintervalTooSmall<Points>( left_half_problem_bin,
                                           right_half_problem_bin ) )
    {
      THROW_EXCEPTION( Utility::IntegratorException,
                       "Error: The minimum subinterval size was reached "
                       "before convergence!" );
    }

    // Check if the roundoff error is too high (store test integrals in counters)
    try{
      this->checkRoundoffError( problem_bin, 
                                left_half_problem_bin, 
                                right_half_problem_bin,    
                                round_off_failed_test_counter_a,
                                round_off_failed_test_counter_b,
                                subinterval+1 );
    }
    EXCEPTION_CATCH_RETHROW( Utility::IntegratorException,
                             "Error: Convergence could not be achieved due to "
                             "roundoff error!" );

    // Improve previous approximations to integral
    initial_bin.getIntegral() += left_half_problem_bin.getIntegral() +
      right_half_problem_bin.getIntegral() - problem_bin.getIntegral();

    // Improve previous approximations to error
    initial_bin.getAbsoluteError() += left_half_problem_bin.getAbsoluteError()+
      right_half_problem_bin.getAbsoluteError() -
      problem_bin.getAbsoluteError();
    
    // Get the error tolerance
    tolerance =
      std::max( d_absolute_error_tol,
                d_relative_error_tol*fabs(initial_bin.getIntegral()) );

    // Check if integral has converged
    if( absolute_error <= tolerance )
      break;
    
    // Subinterval limit has been reached
    else if( subinterval+1 == d_subinterval_limit )
    {
      THROW_EXCEPTION( Utility::IntegratorException,
                       "Error: Maximum number of subdivisions reached before "
                       "convergence!" );
    }
    
    // Continue to next subinterval
    else
    {
      // Add the new subinterval bins to the priority queue
      bin_queue.push( left_half_problem_bin );
      bin_queue.push( right_half_problem_bin );
    }
  }
}

// Bisect and integrate the given bin interval
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<int Points, typename Functor>
void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::bisectAndIntegrateBinInterval( 
                                     Functor& integrand, 
                                     const QuadratureBinType& full_bin,
                                     QuadratureBinType& left_half_bin,
                                     QuadratureBinType& right_half_bin ) const
{
  // Bisect the full bin 
  ArgQuantity midpoint =
    (full_bin.getLowerLimit() + full_bin.getUpperLimit())/2;

  // Set the left half bin
  left_half_bin = QuadratureBinType( full_bin.getLowerLimit(), midpoint );
  
  // Set the right half bin
  right_half_bin = QuadratureBinType( midpoint, full_bin.getUpperLimit() );

  // Integrate over the left half bin
  this->integrateWithPointRule<Points>( integrand, left_half_bin );

  // Integrate over the right half bin
  this->integrateWithPointRule<Points>( integrand, right_half_bin );
};

// Test if subinterval is too small
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<int Points>
inline bool GaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::subintervalTooSmall( 
                                     const QuadratureBinType& left_bin,
                                     const QuadratureBinType& right_bin ) const
{
  // Calculate the epsilon value used for testing valid bin sizes
  unsigned c = Points/10;
  ArgQuantity epsilon = static_cast<FloatType>(1000*c)*ArgQT::epsilon();

  // Get max and min values for testing valid bin sizes
  ArgQuantity max = getMax( fabs( left_bin.getLowerLimit() ),
                            fabs( right_bin.getUpperLimit() ) );
  ArgQuantity min = static_cast<FloatType>(10000)*ArgQT::min();

  if ( max <= (ArgQT::one()+epsilon)*(fabs( right_bin.getLowerLimit() )+min) )
    return true;
  else
    return false;
};

// Check the roundoff error
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::checkRoundoffError( 
                                       const QuadratureBinType& full_bin, 
                                       const QuadratureBinType& left_half_bin, 
                                       const QuadratureBinType& right_half_bin,
                                       int& round_off_failed_test_counter_a,
                                       int& round_off_failed_test_counter_b,
                                       const int number_of_iterations ) const
{
  if( left_half_bin.getIntegralAsc() != left_half_bin.getAbsoluteError() &&
      right_half_bin.getIntegralAsc() != right_half_bin.getAbsoluteError() )
  {
    IntegralQuantity refined_integral = left_half_bin.getIntegral() +
      right_half_bin.getIntegral();
    
    IntegralQuantity refined_error = left_half_bin.getAbsoluteError() +
      right_half_bin.getAbsoluteError();

    IntegralQuantity integral_delta =
      fabs(full_bin.getIntegral() - refined_integral);

    // Check for a diverging or slowly converging integral and error
    if( integral_delta <= fabs(refined_integral)/100000 && 
        refined_error >= 99*full_bin.getAbsoluteError()/100 )
      ++round_off_failed_test_counter_a;
       
    // Check for a diverging error estimate
    if( number_of_iterations >= 10 && error_12 > bin.error )
      ++round_off_failed_test_counter_b;
  }

  // Check if the counters have reached error limits
  TEST_FOR_EXCEPTION( round_off_failed_test_counter_a >= 6,
                      Utility::IntegratorException,
                      "Error: Roundoff error prevented tolerance from being "
                      "achieved (iteration " << number_of_iterations << ")!" );

  TEST_FOR_EXCEPTION( round_off_failed_test_counter_b >= 20,
                      Utility::IntegratorException,
                      "Error: Roundoff error prevented tolerance from being "
                      "achieved (iteration " << number_of_iterations << ")!" );
                      
};

// Conduct the first iteration of the adapative integration
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<typename Functor, typename ArrayType>
bool UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::initializeBinsWynnEpsilon(
                                 Functor& integrand,
                                 const ArrayType& points_of_interest,
                                 BinArray& bin_array,
                                 IntegralQuantity& integral,
                                 IntegralQuantity& absolute_error,
                                 IntegralQuantity& total_absolute_error,
                                 IntegralQuantity& tolerance,
                                 int& ksgn ) const
{
  // Initialize the bin array
  bin_array.resize( d_subinterval_limit );

  // Initialize the bin order array
  unsigned number_of_bins = points_of_interest.size() - 1;
  
  // Keep track of which bins need their errors rescaled
  std::vector<bool> rescale_bin_error( number_of_bins );

  // Initialize the integral
  integral = IntegralQT::zero();

  // Initialize the integral abs
  IntegralQuantity integral_abs = IntegralQT::zero();
  
  // Initialize the absolute error
  absolute_error = IntegralQT::zero();

  // Initialize the total absolute error
  total_absolute_error = IntegralQT::zero();

  // Compute the integral between the points of interest
  for( int i = 0; i < number_of_bins; i++ )
  {
    bin_array[i] = QuadratureBinType( points_of_interest[i],
                                      points_of_interest[i+1] );

    this->integrateWithPointRule<21>( integrand, bin_array[i] );

    bin_array[i].setLevel( 0 );

    // Update the integral, integral abs, and absolute error values
    integral += bin_array[i].getIntegral();
    integral_abs += integral_abs;
    absolute_error += bin_array[i].getAbsoluteError();

    // Check if the bin error needs to be rescaled
    if( bin_array[i].getAbsoluteError() == bin_array[i].getIntegralAsc() &&
        bin_array[i].getAbsoluteError() != IntegralQT::zero() )
      rescale_bin_error[i] = true;
    else
      rescale_bin_error[i] = false;
  }

  // Calculate the convergence parameters
  tolerance = std::max( d_absolute_error_tol,
                        d_relative_error_tol*fabs(integral) );

  IntegralQuanttity round_off =
    100*std::numeric_limits<FloatType>::epsilon()*integral_abs;

  // Check for convergence
  if( absolute_error <= tolerance )
    return true;

  // Check for bad behavior
  else if( absolute_error <= round_off && absolute_error > tolerance )
  {
    THROW_EXCEPTION( Utility::IntegratorException,
                     "Error: Cannot reach tolerance because of roundoff "
                     "error during the first iteration!" );
  }
  
  // Normal behavior - not converged yet
  else
  {
    // Rescale the error approximations for integrals between the points of
    // interest
    for( int i = 0; i < number_of_bins; i++ )
    {
      if( rescale_bin_error[i] ) 
        bin_array[i].setAbsoluteError( absolute_error );

      total_absolute_error += bin_array[i].getAbsoluteError();
    }

    // Set the ksgn value (used for convergence testing)
    if( fabs(integral) >= IntegralQT::one() - 50*round_off )
      ksgn = 1;
    else
      ksgn = -1;
      
    return false;
  }
}

// Conduct the other iterations of the adaptive Wynn Epsilon integration
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<typename Functor>
bool UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::integrateAdaptivelyWynnEpsilonIterate(
                                      Functor& integrand,
                                      BinArray& bin_array,
                                      IntegralQuantity& integral,
                                      IntegralQuantity& absolute_error,
                                      IntegralQuantity& total_absolute_error,
                                      const IntegralQuantit& initial_tolerance,
                                      const unsigned number_of_initial_bins,
                                      const int ksgn ) const
{
  // Make sure the bin array has been initialized properly
  testPrecondition( bin_array.size() == d_subinterval_limit );

  // Reset the absolute error
  absolute_error = IntegralQT::max();
  
  // Set the number of bins
  unsigned number_of_bins = number_of_initial_bins; 

  // Sort bins from greatest to lowest error (don't include the bins that
  // have not been initialized yet)
  std::sort( bin_array.rend()-number_of_bins, bin_array.rend() );
    
  // Initialize the bin order array (keeps the errors sorted)
  std::vector<int> bin_order( number_of_bins );

  for( int i = 0; i < number_of_bins; ++i )
    bin_order[i] = i;

  // This index keeps track of the current problem bin (used in conjunction
  // with the bin_order array
  int nr_max = 0;

  // The maximum number of bisections that may be done on a subinterval
  int max_level = 1;

  // Determines if integrand behavior is currently allowing extrapolation
  bool extrapolation_allowed = true;

  // Determines when at least one bin has reach the max number of bisections
  // and requires an extrapolation.
  bool extrapolation_required = false;

  // Initialize the extrapolated bin integral array
  std::vector<FloatType> bin_extrapolated_integral( 52 );

  // Integral is set during the initialization step
  bin_extrapolated_integral[0] = integral;

  // Keep track of the last three extrapolated integrals
  std::vector<FloatType> last_three_integrals( 3 );
    
  // Keeps track of the number of extrapolations that have been done
  int number_of_extrapolated_calls = 0;

  // Keeps track of the number of extrapolated intervals
  int number_of_extrapolated_intervals = 0;

  // Keep track of th extrapolated tolerance
  IntegralQuantity extrapolated_tolerance = initial_tolerance;
  
  int ierro = 0;
  
  T error_correction = (T)0;

  // Keep track of the error over the "large" bins - all bins that have
  // a level < max_level
  IntegralQuantity error_over_large_bins = total_absolute_error;
  
  // These test counters will be used to determine if the desired tolerance
  // can be acheived
  int round_off_failed_test_counter_a = 0;
  int round_off_failed_test_counter_b = 0;
  int round_off_failed_test_counter_c = 0;
  int ktmin = 0;
  
  for( ; number_of_bins < d_subinterval_limit; ++number_of_bins )
  {
    // Get the problem bin (bin with the largest error)
    // Note: The bins were sorted from highest error to lowest error
    //       during the initialization step.
    ExtrapolatedQuadratureBinType problem_bin = bin_array[bin_order[nr_max]];
    
    ExtrapolatedQuadratureBinType
      left_half_prolem_bin, right_half_problem_bin;
    
    // Always use the 21 point rule to integrate the two bin halfs
    this->bisectAndIntegrateBinInterval<21>( integrand,
                                             problem_bin,
                                             left_half_problem_bin,
                                             right_half_problem_bin );
    
    left_half_problem_bin.setLevel( problem_bin.getLevel() + 1 );
    right_half_problem_bin.setLevel( problem_bin.getLevel() + 1 );
    
    // Check if the subinterval has gotten too small
    if( this->subintervalTooSmall<Points>( left_half_problem_bin,
                                           right_half_problem_bin ) )
    {
      THROW_EXCEPTION( Utility::IntegratorException,
                       "Error: The minimum subinterval size was reached "
                       "before convergence!" );
    }
    
    // Check that the roundoff error is not too high
    try{
      this->checkRoundoffError( problem_bin, 
                                left_half_problem_bin, 
                                right_half_problem_bin,    
                                round_off_failed_test_counter_a,
                                round_off_failed_test_counter_b,
                                round_off_failed_test_counter_c,
                                extrapolation_required,
                                number_of_bins );
    }
    EXCEPTION_CATCH_RETHROW( Utility::IntegratorException,
                             "Error: Convergence could not be achieved due "
                             "to roundoff error!" );

    // Improve previous approximations to the integral and error 
    integral += left_half_problem_bin.getIntegral() +
      right_half_problem_bin.getIntegral() - problem_bin.getIntegral();
    
    total_absolute_error += left_half_problem_bin.getAbsoluteError() +
      right_half_problem_bin.getAbsoluteError() -
      problem_bin.getAbsoluteError();
    

    // Update and sort bin order
    this->sortBins( left_half_problem_bin,
                    right_half_problem_bin,
                    bin_order, 
                    bin_array,
                    nr_max,
                    number_of_bins );
    
    // Get the error tolerance
    IntegralQuantity tolerance = 
      std::max( d_absolute_error_tol,
                d_relative_error_tol*fabs(initial_bin.getIntegral()) );

    // Check if the integral has converged
    if( total_absolute_error <= tolerance )
      break;

    // Subinterval limit has been reached
    else if( number_of_bins+1 == d_subinterval_limit )
    {
      THROW_EXCEPTION( Utility::IntegratorException,
                       "Error: Maximum number of subdivisions reached before "
                       "convergence!" );
    }

    // Extrapolation is currently allowed 
    else if( extrapolation_allowed )
    {
      // Update the error over "large" bins - bins that have not
      // reach the max subdivision limit yet.
      error_over_large_bins -= problem_bin.getAbsoluteError();

      if( left_half_problem_bin.getLevel() + 1 <= max_level )
      {
        error_over_large_bins += left_half_problem_bin.getAbsoluteError()+
          right_half_problem_bin.getAbsoluteError();
      }

      // Check if an extrapolation is required (at least one bin has
      // reach the max subdivision level)
      extrapolation_required =
        this->checkIfExtrapolationIsRequired( bin_array,
                                              bin_order,
                                              extrapolation_required,
                                              max_level,
                                              nr_max );

      // Check if there are any bins that have not reached the max
      // subdivision limit - extrapolation can not be done until there
      // are no more "large" bins.
      if( extrapolation_required )
      {
        const bool all_bins_ready =
          this->checkIfAllBinsReadyForExtrapolation( bin_array,
                                                     bin_order,
                                                     error_over_large_bins,
                                                     extrapolated_tolerance,
                                                     number_of_bins,
                                                     max_level,
                                                     nr_max );

        // Conduct the extrapolation
        if( all_bins_ready )
        {
          ++number_of_extrapolated_intervals;
          
          bin_extrapolated_integral[number_of_extrapolated_intervals] =
            integral;

          // Increate the max bisection level and restart
          if( number_of_extrapolated_intervals < 2 )
          {
            nr_max = 0;
            extrapolation_required = false;
            max_level++;
            error_over_large_bins = total_absolute_error;
            continue;
        }

        IntegralQuantity extrapolated_integral = IntegralQT::zero();
        IntegralQuantity extrapolated_error = IntegralQT::zero();
          
        this->getWynnEpsilonAlgorithmExtrapolation( 
                                              bin_extrapolated_integral,
                                              last_three_integrals,
                                              extrapolated_integral, 
                                              extrapolated_error,  
                                              number_of_extrapolated_intervals,
                                              number_of_extrapolated_calls );

        ktmin++;
        
        TEST_FOR_EXCEPTION( ktmin > 5 &&
                            absolute_error < total_absolute_error/1000, 
                            Utility::IntegratorException,
                            "Error: The integral is probably divergent "
                            "(or slowly convergent)!" );
        
        if( extrapolated_error < absolute_error )
        {
          ktmin = 0;
          absolute_error = extrapolated_error;
          integral = extrapolated_integral;
          error_correction = error_over_large_bins;
          extrapolated_tolerance = 
            this->getMax( d_absolute_error_tol, 
                          d_relative_error_tol*fabs( extrapolated_integral ) );

          if( absolute_error <= extrapolated_tolerance )
            break;
        }
        
        // Prepare bisection of the smallest interval.
        if( number_of_extrapolated_intervals == 0 )
          no_extrapolation_allowed = true;

        nr_max = 0;
        extrapolation_required = false;
        max_level++;
        error_over_large_bins = total_error;
      }
    }
  } // end main for loop

  //  Set final integral and error estimate.
  if( absolute_error == std::numeric_limits<T>::max() )
  {
    //  Compute global integral sum.  
    T long_integral = (T)0;
    std::vector<int>::reverse_iterator j =  bin_order.rbegin();
    // Sum integral over all bins
    for( j; j != bin_order.rend(); j++ )
    {
      bin = bin_array[*j];
      long_integral += bin.integral;
    }
    integral = long_integral;
    absolute_error = total_error;

    return;
  }

  // Test on divergence.
  if( ksgn == (-1) && this->getMax( fabs(integral), fabs(total_area) ) <=
      (total_area_abs/(T)100) )
  {
    return;
  }

  TEST_FOR_EXCEPTION( (1/(T)100) > integral/total_area ||
                      integral/total_area > (T)100 ||
                      total_error > fabs( total_area ), 
                      Utility::IntegratorException,
                      "Error: The input is invalid! Because the absolute "
                      "error tol < 0 and the relative error tol < 0, the "
                      "integral and absolute_error are set to zero." );

  return;
}

// Sort the bin order from highest to lowest error 
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::sortBins(
                           const ExtrapolatedQuadratureBinType& left_half_bin,
                           const ExtrapolatedQuadratureBinType& right_half_bin,
                           std::vector<int>& bin_order,
                           BinArray& bin_array,
                           int& nr_max,
                           const int number_of_bins ) const
{
  // Make sure the bin order array is valid
  testPrecondition( bin_order.size() == number_of_bins );

  // Record the larger and smaller errors of the two new bins
  IntegralQuantity larger_new_bin_error;
  IntegralQuantity smaller_new_bin_error;

  // Get the bin index corresponding to the bin with the largest error
  int bin_with_larger_error = bin_order[nr_max];

  // Add the new bins to the bin array
  if ( left_half_bin.getAbsoluteError() <= right_half_bin.getAbsoluteError() )
  {
    bin_array[bin_with_larger_error] = right_half_bin;
    bin_array[number_of_bins] = left_half_bin;

    larger_new_bin_error = right_half_bin.getAbsoluteError();
    smaller_new_bin_error = left_half_bin.getAbsoluteError();
  }
  else
  {
    bin_array[bin_with_larger_error] = left_half_bin;
    bin_array[number_of_bins] = right_half_bin;

    larger_new_bin_error = left_half_bin.getAbsoluteError();
    smaller_new_bin_error = right_half_bin.getAbsoluteError();
  }

  // Remove the old bin index from the bin order array 
  bin_order.remove( nr_max ); 

  // This part of the routine is only executed if, due to a
  // difficult integrand, subdivision increased the error
  // estimate. in the normal case the insert procedure should
  // start after the nr_max-th largest error estimate.
  int original_nr_max = nr_max;
  
  while( nr_max > 0 && larger_new_bin_error >
         bin_array[bin_order[nr_max-1]].getAbsoluteError() )
  {
    // Reduce nr_max if the bin above it has larger error
    --nr_max; 
  }

  // Start insert of the bin_with_larger_error at the reduced nr_max
  int start_bin_index;
  
  if( original_nr_max > nr_max )
    start_bin_index = nr_max;
  
  // Start insert of the bin_with_larger_error right after the
  // bin_with_larger_error
  else
    start_bin_index = original_nr_max;

  // Compute the number of elements in the list to be maintained
  // in descending order. This number depends on the number of
  // subdivisions still allowed.
  std::vector<int>::iterator max_bin_it;
  
  if( (d_subinterval_limit/2+2) < bin_order.size()-1 )
    max_bin_it = bin_order.begin() + (d_subinterval_limit - bin_order.size());
  else
    max_bin_it = bin_order.end();

  // Insert the bin with the larger error by traversing the list top-down.
  std::vector<int>::iterator large_bin_it = bin_order.begin()+start_bin_index;
  
  while( large_bin_it != max_bin_it &&
         larger_new_bin_error < bin_array[*large_bin_it].getAbsoluteError() )
  {
    ++large_bin_it;
  }

  bin_order.insert( large_bin_it, bin_with_larger_error );

  //  Insert smaller_bin_error by traversing the list bottom-up.
  std::vector<int>::iterator small_bin_it = max_bin_it;
  
  while( small_bin_it != large_bin_it && 
         bin_array[bin_order[*small_bin_it]].getAbsoluteError() <
         smaller_new_bin_error )
  {
    --small_bin_it;
  }
  
  bin_order.insert( small_bin_it+1, number_of_bins );
}

// Check if an extrapolation is required on the current worst bin
// Note: If extrapolation_required is already TRUE, this method does not
//       need to check anything. If extraplation_required is FALSE, the
//       current worst bin will be checked to see if it has reached the
//       bisection limit. If it has, TRUE will be returned and nr_max will
//       be reset so that a new problem bin in need of more bisections can
//       be found.
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::checkIfExtrapolationIsRequired(
                                       const BinArray& bin_array,
                                       const std::vector<int>& bin_order_array,
                                       const bool extrapolation_required,
                                       const int max_level,
                                       int& nr_max ) const
{
  // Extrapolation is already required in at least one bin
  if( extrapolation_required )
    return true;

  // Extrapolation is not required yet - check if the current worst bin
  // requires extrapolation
  else
  {
    const ExtrapolatedQuadratureBinType& problem_bin =
      bin_array[bin_order_array[nr_max]];
    
    // The bin must be bisected further before extrapolating
    if( problem_bin.getLevel() + 1 <= max_level )
      return false;
    
    // The bin has been bisected to the maximum level - extrapolation required
    else
    {
      // We know that the bin at bin_order_array[nr_max] requires
      // extrapolation. Reset nr_max so that we can search for other bins
      // that require extrapolation (later).
      nr_max = 1;

      return true;
    }
  }
}

// Check if all bins are ready for extrapolation
// Note: If this returns false, nr_max will also be set to a bin that is
//       still in need of more bisections.
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::checkIfAllBinsReadyForExtrapolation(
                                     const BinArray& bin_array,
                                     const std::vector<int>& bin_order_array,
                                     const IntegralType& error_over_large_bins,
                                     const IntegralType& tolerance,
                                     const unsigned number_of_bins,
                                     const int max_level,
                                     int& nr_max ) const
{
  bool all_bins_ready = true;

  // Check if the large bin error is still large
  if( error_over_large_bins > tolerance )
  {
    int size = number_of_bins;
          
    if( number_of_bins >  2 + d_subinterval_limit/2 )
      size = d_subinterval_limit + 3 - number_of_bins;

    bool still_have_large_bins = false;

    // Check if there are any bins that are still large
    for( int k = nr_max; k < size; k++ )
    {
      const ExtrapolatedQuadratureBinType& bin = bin_array[bin_order[nr_max]];
      
      if( bin.getLevel() + 1 <= max_level )
      {
        still_have_large_bins = true;

        break;
      }
      else
        ++nr_max;    
    }    
          
    if( still_have_large_bins )
      all_bins_ready = false;
  }

  return all_bins_ready;
}

// get the Wynn Epsilon-Algoirithm extrapolated value
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::extrapolateWithWynnEpsilonAlgorithm( 
                      std::vector<IntegralQuantity>& bin_extrapolated_integral,
                      std::vector<IntegralQuantity>& last_three_integrals, 
                      IntegralQuantity& extrapolated_integral, 
                      IntegralQuantity& extrapolated_error,  
                      int& number_of_extrapolated_intervals,
                      int& number_of_extrapolated_calls  ) const
{
  // Make sure the bin extrapolated integral array is valid
  testPrecondition( bin_extrapolated_integral.size() == 52 );
  // Make sure the last three integrals array is valid
  testPrecondition( last_three_integrals.size() == 3 );
  // Make sure the number of extrapolated intervals is valid
  testPrecondition( number_of_extrapolated_intervals > 1 );
  // Make sure the number of extrapolated calls is valid
  testPrecondition( number_of_extrapolated_calls >= 0 );
 
  // update the number of extrapolated calls
  ++number_of_extrapolated_calls;

  extrapolated_error = std::numeric_limits<T>::max();
  extrapolated_integral = bin_extrapolated_integral[number_of_extrapolated_intervals];

  if ( number_of_extrapolated_intervals < 2 )
  {
    extrapolated_error = getMax( extrapolated_error, 
        std::numeric_limits<T>::epsilon()*fabs(extrapolated_integral)/(T)2 );
    return;
  }

  int extrapolated_interval_limit = 50;

  bin_extrapolated_integral[number_of_extrapolated_intervals+2] = 
    bin_extrapolated_integral[number_of_extrapolated_intervals];

  int new_element = number_of_extrapolated_intervals/2;

  bin_extrapolated_integral[number_of_extrapolated_intervals] = 
     std::numeric_limits<T>::max();

  int original_number = number_of_extrapolated_intervals;
  int k1 = number_of_extrapolated_intervals;

  for ( int i = 0; i < new_element; i++ )
  {
    int k2 = k1-1;
    int k3 = k1-2;

    T integral = bin_extrapolated_integral[k1+2];
    T e0 = bin_extrapolated_integral[k3];
    T e1 = bin_extrapolated_integral[k2];
    T e2 = integral;

    // Get error and tolerance estimate between e2 and e1
    T delta2 = e2 - e1;
    T error2 = fabs(delta2);
    T tolerance2 = getMax( fabs(e2), fabs(e1) )*
        std::numeric_limits<T>::epsilon();

    // Get error and tolerance estimate between e1 and e0
    T delta3 = e1 - e0;
    T error3 = fabs(delta3);
    T tolerance3 = getMax( fabs(e1), fabs(e0) )*
        std::numeric_limits<T>::epsilon();

    // If e0, e1 and e2 are equal to within machine accuracy, convergence is assumed.
    if ( error2 <= tolerance2 && error3 <= tolerance3 )
    {
      extrapolated_integral = integral;
      extrapolated_error = error2+error3;
      extrapolated_error = getMax( extrapolated_error, 
          (std::numeric_limits<T>::epsilon()*fabs(extrapolated_error))/(T)2 );
     return;
    }

    T e3 = bin_extrapolated_integral[k1];
    bin_extrapolated_integral[k1] = e1;

    // Get error and tolerance estimate between e1 and e3
    T delta1 = e1 - e3;
    T error1 = fabs(delta1);
    T tolerance1 = getMax( fabs(e1), fabs(e3) )*
            std::numeric_limits<T>::epsilon();
    
    /* If two elements are very close to each other, omit a part of the table 
     * by adjusting the value of number_of_extrapolated_intervals.
     */
    if ( error1 <= tolerance1 || error2 <= tolerance2 || error3 <= tolerance3 )
    {
      number_of_extrapolated_intervals = 2*i;
      break;
    }

    T ss = (1/delta1 + 1/delta2) - 1/delta3;

    /* Test to detect irregular behavior in the table, and eventually omit a
     * part of the table adjusting the value of number_of_extrapolated_intervals.
     */
    if ( fabs ( ss*e1 ) <= 1/(T)10000 ) 
    {
      number_of_extrapolated_intervals = 2*i;
      break;
    }

    /* Compute a new element and eventually adjust the value of
     * integral. 
     */
    integral = e1 + 1/ss;
    bin_extrapolated_integral[k1] = integral;
    k1 -= 2;

    T error = error2 + fabs(integral - e2) + error3;

    if ( error <= extrapolated_error )
    {
      extrapolated_error = error;
      extrapolated_integral = integral;
    }
  }

  // Shift the table

  if ( number_of_extrapolated_intervals == extrapolated_interval_limit )
    number_of_extrapolated_intervals = 2*(extrapolated_interval_limit/2);

  int ib;

  if ( original_number % 2 == 1 )
  {
    ib = 1;
  }
  else
  {
    ib = 0;
  }

  int ie = new_element + 1;

  for ( int i = 0; i < ie; i++ )
  {
    bin_extrapolated_integral[ib] = bin_extrapolated_integral[ib+2];
    ib += 2;
  }

  if ( original_number != number_of_extrapolated_intervals ) 
  {
    for ( int i = 0; i < number_of_extrapolated_intervals; i++ )
    {
      bin_extrapolated_integral[i]= 
        bin_extrapolated_integral[original_number - number_of_extrapolated_intervals + i];
    }
  }
  
  if ( number_of_extrapolated_calls < 4 )
  {
    last_three_integrals[number_of_extrapolated_calls-1] = extrapolated_integral;
    extrapolated_error = std::numeric_limits<T>::max();
  }
  else
  {
    extrapolated_error = 
      fabs( extrapolated_integral - last_three_integrals[2] ) +
      fabs( extrapolated_integral - last_three_integrals[1] ) +
      fabs( extrapolated_integral - last_three_integrals[0] );

    last_three_integrals[0] = last_three_integrals[1];
    last_three_integrals[1] = last_three_integrals[2];
    last_three_integrals[2] = extrapolated_integral;
  }

  extrapolated_error = std::max ( 
    extrapolated_error,
    (std::numeric_limits<T>::epsilon()*fabs( extrapolated_integral ))/(T)2 );

  return;
};  

// check the roundoff error
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
void UnitAwareGaussKronrodIntegrator<FloatType,ArgUnit,IntegrandUnit>::checkRoundoffError( 
                           const ExtrapolatedQuadratureBinType& full_bin, 
                           const ExtrapolatedQuadratureBinType& left_half_bin, 
                           const ExtrapolatedQuadratureBinType& right_half_bin,
                           int& round_off_failed_test_counter_a,
                           int& round_off_failed_test_counter_b,
                           int& round_off_failed_test_counter_c,
                           const bool extrapolate, 
                           const int number_of_iterations ) const
{
  if( left_half_bin.getIntegralAsc() != left_half_bin.getAbsoluteError() &&
      right_half_bin.getIntegralAsc() != right_half_bin.getAbsoluteError() )
  {
    
    IntegralQuantity refined_integral =
      bin_1.getIntegral() + bin_2.getIntegral();
       
    IntegralQuantity refined_error =
      bin_1.getAbsoluteError() + bin_2.getAbsoluteError();

    IntegralQuantity integral_delta =
      fabs(full_bin.getIntegral() - refined_integral);

    // Check for a diverging or slowly converging integral and error
    if( integral_delta <= fabs(refined_integral)/100000 && 
        refined_error >= 99*bin.getAbsoluteError()/100 )
    {
      if( extrapolate ) 
        ++round_off_failed_test_counter_b;
      else
        ++round_off_failed_test_counter_a; 
    }

    // Check for a diverging error estimate
    if( number_of_iterations >= 10 && refined_error > bin.getAbsoluteError() )
      ++round_off_failed_test_counter_c;
  }
  
  // Check if the counters have reached error limits
  TEST_FOR_EXCEPTION( round_off_failed_test_counter_a +
                      round_off_failed_test_coutner_b >= 10,
                      Utility::IntegratorException,
                      "Error: Roundoff error prevented tolerance from being "
                      "achieved (iteration " << number_of_iterations << ")!" );
  
  TEST_FOR_EXCEPTION( round_off_failed_test_counter_b >= 5,
                      Utility::IntegratorException,
                      "Error: Extremely bad integrand behavior occurs at some "
                      "points of the integration interval! The tolerance "
                      "cannot be achieved (iteration "
                      << number_of_iterations << ")!" );

  TEST_FOR_EXCEPTION( round_off_failed_test_counter_c >= 20,
                      Utility::IntegratorException,
                      "Error: Roundoff error prevented tolerance from being "
                      "achieved (iteration " << number_of_iterations << ")!" );
};

} // end Utility namespace

#endif // end UTILITY_GAUSS_KRONROD_INTEGRATOR_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_GaussKronrodIntegrator_def.hpp
//---------------------------------------------------------------------------//
