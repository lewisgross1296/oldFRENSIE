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

// Integrate the function adaptively with BinQueue
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
                                       ArgQuantity lower_limit, 
                                       ArgQuantity upper_limit,
                                       IntegralQuantity& result,
				       IntegralQuantity& absolute_error ) const
{
  // Make sure the integration limits are valid
  testPrecondition( lower_limit <= upper_limit );
  
  // Create the bin that encompasses the entire integration range
  QuadratureBinType bin( lower_limit, upper_limit );
  
  bool converged;

  try{
    converged = this->integrateAdaptivelyInitialIteration( integrand, bin );
  }
  EXCEPTION_CATCH_RETHROW( Utility::IntegratorException,
                           "Error: The initial iteration of the adaptive "
                           "integration failed!" );
  
  // Check for fast convergence
  if( converged )
  {
    result = bin.getResult();
    absolute_error = bin.getAbsoluteError();
    
    return;
  }
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
      this->integrateAdaptivelyIterate( integrand,
                                        bin,
                                        result,
                                        absolute_error );
    }
    EXCEPTION_CATCH_RETHROW( Utility::IntegratorException,
                             "Error: The integration failed during the "
                             "subinterval iterations."
  }
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
                                           ArgQuantity lower_limit, 
                                           ArgQuantity upper_limit,
                                           IntegralQuantity& result,
                                           IntegralQuantity& absolute_error,
                                           IntegralQuantity& result_abs, 
                                           IntegralQuantity& result_asc ) const
{
  // Make sure the integration limits are valid
  testPrecondition( lower_limit <= upper_limit );

  // Get the quadrature set traits for the requested point rule
  typedef GaussKronrodQuadratureSetTraits<ArgUnit,FloatType> GKQST;
  
  if( lower_limit < upper_limit )
  {
    // Calculate the midpoint
    ArgQuantity midpoint = (upper_limit + lower_limit)/2;

    // Calculate the half the length between the upper and lower integration
    // limits
    ArgQuantity half_length = (upper_limit - lower_limit )/(T)2;
    ArgQuantity abs_half_length = fabs( half_length );

    // Get number of Kronrod weights
    int number_of_weights = GKQST::getKronrodWeights().size();

    Teuchos::Array<IntegrandQuantity>
      integrand_values_lower( number_of_weights );
    Teuchos::Array<IntegrandQuantity>
      integrand_values_upper( number_of_weights );
    Teuchos::Array<IntegrandQuantity>
      integrand_values_sum( number_of_weights );

    // Estimate Kronrod and absolute value integral for all but last weight
    IntegralQuantity kronrod_result = IntegralQT::zero();
    result_abs = IntegralQT::zero();
    
    for( int j = 0; j < number_of_weights-1; j++ )
    {  
      this->calculateQuadratureIntegrandValuesAtAbscissa( 
                                               integrand,
                                               GKQST::getKronrodAbscissae()[j],
                                               half_length,
                                               midpoint,
                                               integrand_values_lower[j],
                                               integrand_values_upper[j] );

        
      integrand_values_sum[j] = 
        integrand_values_lower[j] + integrand_values_upper[j];
      
      kronrod_result +=
        GKQST::getKronrodWeights()[j]*integrand_values_sum[j];
      
      result_abs += GKQST::getKronrodWeights()[j]*
        (fabs( integrand_values_lower[j] ) +
         fabs( integrand_values_upper[j] ) );
    };

    // Integrand at the midpoint
    IntegrandQuantity integrand_midpoint = integrand( midpoint );

    // Estimate Kronrod integral for the last weight
    IntegralQuantity kronrod_result_last_weight = integrand_midpoint*
      GKQST::getKronrodWeights().back();

    // Update Kronrod estimate and absolute value with last weight
    kronrod_result += kronrod_result_last_weight;
    result_abs += fabs( kronrod_result_last_weight );

    // Calculate final integral result and absolute value
    result = kronrod_result*half_length;
    result_abs *= abs_half_length;

    // Calculate the mean kronrod result
    T mean_kronrod_result = kronrod_result/(T)2;

    // Estimate the result asc for all but the last weight
    result_asc = (T)0;
    for ( int j = 0; j < number_of_weights - 1; j++ )
      {  
        result_asc += GaussKronrodQuadratureSetTraits<Points>::kronrod_weights[j]*
          ( fabs( integrand_values_lower[j] - mean_kronrod_result ) +
            fabs( integrand_values_upper[j] - mean_kronrod_result ) );
      };

    // Estimate the result asc for the last weight
    result_asc += GaussKronrodQuadratureSetTraits<Points>::kronrod_weights[number_of_weights-1]*
        fabs( integrand_midpoint - mean_kronrod_result );

    // Calculate final result acx
    result_asc *= abs_half_length;

    // Estimate Gauss integral
    T gauss_result = (T)0;

    for ( int j = 0; j < (number_of_weights-1)/2; j++ )
      {
        int jj = j*2 + 1;
        gauss_result += integrand_values_sum[jj]*
            GaussKronrodQuadratureSetTraits<Points>::gauss_weights[j];
      };

    // Update Gauss estimate with last weight if needed
    if ( number_of_weights % 2 == 0 )
    {
      gauss_result += integrand_midpoint*
        GaussKronrodQuadratureSetTraits<Points>::gauss_weights[number_of_weights/2 - 1];
    }

    // Estimate error in integral 
    absolute_error = fabs( ( kronrod_result - gauss_result ) * half_length );
    rescaleAbsoluteError( absolute_error, result_abs, result_asc);

  }
  else if( lower_limit == upper_limit )
  {
    result = IntegralQT::zero();
    absolute_error = IntegralQT::zero();
  }
  else // invalid limits
  {
    THROW_EXCEPTION( Utility::IntegratorException,
		     "Invalid integration limits: " << lower_limit << " !< "
		     << upper_limit << "." );
  }
}

// Rescale absolute error from integration
template<typename T>
void UnitAwareGaussKronrodIntegrator<T>::rescaleAbsoluteError( 
    T& absolute_error, 
    T result_abs, 
    T result_asc ) const
{
  if ( result_asc != 0 && absolute_error != 0 )
    {
      T scale = (T)200 * absolute_error/result_asc;

      if ( scale < (T)1 )
      {
        absolute_error = result_asc * pow( scale, 3/(T)2 );
      }  
      else
      {    
        absolute_error = result_asc;
      }  
    };

  if ( result_abs > std::numeric_limits<T>::min() / 
        ( (T)50 * std::numeric_limits<T>::epsilon() ) )
    {
      T min_error = (T)50*std::numeric_limits<T>::epsilon() * result_abs;

      if ( min_error > absolute_error ) 
        {
          absolute_error = min_error;
        }
    };
};

// Sort the bin order from highest to lowest error 
//! \details The error list will be correctly sorted except bin_1 and bin_2
template<typename T>
void UnitAwareGaussKronrodIntegrator<T>::sortBins( 
        Teuchos::Array<int>& bin_order,
        BinArray& bin_array, 
        const ExtrapolatedQuadratureBin<T>& bin_1,
        const ExtrapolatedQuadratureBin<T>& bin_2,
        const int& number_of_intervals,
        int& nr_max ) const
{
  testPrecondition( bin_order.size() == number_of_intervals );

  T larger_error;
  T smaller_error;
  
  int bin_with_larger_error = bin_order[nr_max];

  // append new intervals to bin_array
  if ( bin_1.getAbsoluteError() <= bin_2.getAbsoluteError() )
  {
    bin_array[bin_with_larger_error] = bin_2;
    bin_array[number_of_intervals] = bin_1;

    larger_error = bin_2.getAbsoluteError();
    smaller_error = bin_1.getAbsoluteError();
  }
  else
  {
    bin_array[bin_with_larger_error] = bin_1;
    bin_array[number_of_intervals] = bin_2;

    larger_error = bin_1.getAbsoluteError();
    smaller_error = bin_2.getAbsoluteError();
  }

  // remove old interval from list
  bin_order.remove( nr_max ); 

  /*
   *  This part of the routine is only executed if, due to a
   *  difficult integrand, subdivision increased the error
   *  estimate. in the normal case the insert procedure should
   *  start after the nr_max-th largest error estimate.
   */
  int original_nr_max = nr_max;
  while ( nr_max > 0 && larger_error >
          bin_array[bin_order[nr_max-1]].getAbsoluteError() )
  {
    --nr_max; //reduce nr_max if the bin above it has larger error
  }

  // Start insert of the bin_with_larger_error at the reduced nr_max
  int start_bin;
  if ( original_nr_max > nr_max )
    start_bin = nr_max;
  
  // Start insert of the bin_with_larger_error right after the
  // bin_with_larger_error
  else
    start_bin = original_nr_max;

  /*
   *  Compute the number of elements in the list to be maintained
   *  in descending order. This number depends on the number of
   *  subdivisions still allowed.
   */
  Teuchos::Array<int>::iterator max_bin;
  if ( (d_subinterval_limit/2+2) < bin_order.size()-1 )
    max_bin = bin_order.begin() + ( d_subinterval_limit - bin_order.size() );
  else
    max_bin = bin_order.end();

  Teuchos::Array<int>::iterator large_bin = bin_order.begin()+start_bin;
  while ( large_bin != max_bin &&
          larger_error < bin_array[*large_bin].getAbsoluteError() )
  {
    large_bin++;
  }

  bin_order.insert( large_bin, bin_with_larger_error );
  max_bin;

  //  Insert smaller_bin_error by traversing the list bottom-up.
  Teuchos::Array<int>::iterator small_bin = max_bin;
  while ( small_bin != large_bin && 
          bin_array[bin_order[*small_bin]].getAbsoluteError() < smaller_error )
  {
    --small_bin;
  }
  
  bin_order.insert( small_bin+1, number_of_intervals );
};


// get the Wynn Epsilon-Algoirithm extrapolated value
template<typename T>
void UnitAwareGaussKronrodIntegrator<T>::getWynnEpsilonAlgorithmExtrapolation( 
        Teuchos::Array<T>& bin_extrapolated_result, 
        Teuchos::Array<T>& last_three_results, 
        T& extrapolated_result, 
        T& extrapolated_error,  
        int& number_of_extrapolated_intervals,
        int& number_of_extrapolated_calls  ) const
{
  testPrecondition( number_of_extrapolated_calls >= 0 );
  testPrecondition( number_of_extrapolated_intervals > 1 );
  testPrecondition( last_three_results.size() == 3 );
  testPrecondition( bin_extrapolated_result.size() == 52 );
 
  // update the number of extrapolated calls
  number_of_extrapolated_calls++;

  extrapolated_error = std::numeric_limits<T>::max();
  extrapolated_result = bin_extrapolated_result[number_of_extrapolated_intervals];

  if ( number_of_extrapolated_intervals < 2 )
  {
    extrapolated_error = getMax( extrapolated_error, 
        std::numeric_limits<T>::epsilon()*fabs(extrapolated_result)/(T)2 );
    return;
  }

  int extrapolated_interval_limit = 50;

  bin_extrapolated_result[number_of_extrapolated_intervals+2] = 
    bin_extrapolated_result[number_of_extrapolated_intervals];

  int new_element = number_of_extrapolated_intervals/2;

  bin_extrapolated_result[number_of_extrapolated_intervals] = 
     std::numeric_limits<T>::max();

  int original_number = number_of_extrapolated_intervals;
  int k1 = number_of_extrapolated_intervals;

  for ( int i = 0; i < new_element; i++ )
  {
    int k2 = k1-1;
    int k3 = k1-2;

    T result = bin_extrapolated_result[k1+2];
    T e0 = bin_extrapolated_result[k3];
    T e1 = bin_extrapolated_result[k2];
    T e2 = result;

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
      extrapolated_result = result;
      extrapolated_error = error2+error3;
      extrapolated_error = getMax( extrapolated_error, 
          (std::numeric_limits<T>::epsilon()*fabs(extrapolated_error))/(T)2 );
     return;
    }

    T e3 = bin_extrapolated_result[k1];
    bin_extrapolated_result[k1] = e1;

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
     * result. 
     */
    result = e1 + 1/ss;
    bin_extrapolated_result[k1] = result;
    k1 -= 2;

    T error = error2 + fabs(result - e2) + error3;

    if ( error <= extrapolated_error )
    {
      extrapolated_error = error;
      extrapolated_result = result;
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
    bin_extrapolated_result[ib] = bin_extrapolated_result[ib+2];
    ib += 2;
  }

  if ( original_number != number_of_extrapolated_intervals ) 
  {
    for ( int i = 0; i < number_of_extrapolated_intervals; i++ )
    {
      bin_extrapolated_result[i]= 
        bin_extrapolated_result[original_number - number_of_extrapolated_intervals + i];
    }
  }
  
  if ( number_of_extrapolated_calls < 4 )
  {
    last_three_results[number_of_extrapolated_calls-1] = extrapolated_result;
    extrapolated_error = std::numeric_limits<T>::max();
  }
  else
  {
    extrapolated_error = 
      fabs( extrapolated_result - last_three_results[2] ) +
      fabs( extrapolated_result - last_three_results[1] ) +
      fabs( extrapolated_result - last_three_results[0] );

    last_three_results[0] = last_three_results[1];
    last_three_results[1] = last_three_results[2];
    last_three_results[2] = extrapolated_result;
  }

  extrapolated_error = std::max ( 
    extrapolated_error,
    (std::numeric_limits<T>::epsilon()*fabs( extrapolated_result ))/(T)2 );

  return;
};  


// check the roundoff error
template<typename T>
void UnitAwareGaussKronrodIntegrator<T>::checkRoundoffError( 
                       const QuadratureBin<T>& bin, 
                       const QuadratureBin<T>& bin_1, 
                       const QuadratureBin<T>& bin_2,    
                       const T& bin_1_asc,
                       const T& bin_2_asc,
                       int& round_off_1,
                       int& round_off_2,
                       const int number_of_iterations ) const
{
    if (bin_1_asc != bin_1.error && bin_2_asc != bin_2.error)
    {
       T area_12 = bin_1.result + bin_2.result;
       T error_12 = bin_1.error + bin_2.error;
       T delta = bin.result - area_12;

       if ( fabs (delta) <= (1/(T)100000) * fabs (area_12) && 
            error_12 >= (99/(T)100) * bin.error )
       {
         round_off_1++;
       }
       if ( number_of_iterations >= 10 && error_12 > bin.error )
       {
          round_off_2++;
       }
     }

    TEST_FOR_EXCEPTION( round_off_1 >= 6 || round_off_2 >= 20, 
                        Utility::IntegratorException,
                        "Roundoff error prevented tolerance from being achieved" );
};


// check the roundoff error
template<typename T>
void UnitAwareGaussKronrodIntegrator<T>::checkRoundoffError( 
                       const ExtrapolatedQuadratureBin<T>& bin, 
                       const ExtrapolatedQuadratureBin<T>& bin_1, 
                       const ExtrapolatedQuadratureBin<T>& bin_2,    
                       const T& bin_1_asc,
                       const T& bin_2_asc,
                       int& round_off_1,
                       int& round_off_2,
                       int& round_off_3,
                       const bool extrapolate, 
                       const int number_of_iterations ) const
{
  if( bin_1_asc != bin_1.getAbsoluteError() &&
      bin_2_asc != bin_2.getAbsoluteError() )
    {
       T area_12 = bin_1.getResult() + bin_2.getResult();
       T error_12 = bin_1.getAbsoluteError() + bin_2.getAbsoluteError();
       T delta = bin.getResult() - area_12;

       if( fabs (delta) <= (1/(T)100000) * fabs (area_12) && 
           error_12 >= (99/(T)100) * bin.getAbsoluteError() )
       {
         if( extrapolate ) 
           round_off_2++;
         else
           round_off_1++; 
       }
       if( number_of_iterations >= 10 &&
           error_12 > bin.getAbsoluteError() )
          round_off_3++;
     }

    TEST_FOR_EXCEPTION( 10 <= round_off_1 + round_off_2 || 20 <= round_off_3, 
                        Utility::IntegratorException,
                        "Roundoff error prevented tolerance from being achieved" );

    TEST_FOR_EXCEPTION( 5 <= round_off_2, 
                        Utility::IntegratorException,
                        "Extremely bad integrand behavior occurs at some points "
                        "of the integration interval" );
};

template<typename T>
template<int Points, typename Functor>
bool UnitAwareGaussKronrodIntegrator<T>::integrateAdaptivelyInitialIteration(
						 Functor& integrand, 
						 QuadratureBin<T>& bin ) const
{
  // Perform an integration using the requested point rule over the interval
  // of interest (specified by the bin)
  T result = (T)0;
  T absolute_error = (T)0;
  T result_abs = (T)0;
  T result_asc = (T)0;

  integrateWithPointRule<Points>(
    integrand, 
    bin.getLowerLimit()
    bin.getUpperLimit(), 
    result,
    absolute_error,
    result_abs, 
    result_asc );

  bin.setResultAndAbsoluteError( result, absolute_error );

  // Test on accuracy
  T tolerance = getMax( d_absolute_error_tol,
                        d_relative_error_tol*fabs( bin.result ) );

  // Need IEEE rounding here to match original quadpack behavior
  T volatile_round_off;
  volatile_round_off = (T)50*std::numeric_limits<T>::epsilon()*result_abs;
  T round_off = volatile_round_off;

  TEST_FOR_EXCEPTION( bin.getAbsoluteError() <= round_off &&
                      bin.getAbsoluteError() > tolerance, 
                      Utility::IntegratorException,
                      "cannot reach tolerance because of roundoff error "
                      "on first attempt" );

  // The integral has converged
  if( (bin.getAbsoluteError() <= tolerance &&
       bin.getAbsoluteError() != result_asc) || 
      bin.getAbsoluteError() == (T)0 ) 
    return true;
  // The integral has not converged - more iterations will be required
  else
    return false;
}

// Conduct the other iterations of the adaptive integration
template<int Points, typename Functor>
void UnitAwareGaussKronrodIntegrator<T>::integrateAdaptivelyIterate(
                                                 Functor& integrand,
                                                 QuadratureBin<T>& initial_bin,
                                                 T& result,
                                                 T& absolute_error ) const
{
  BinQueue bin_queue;
  bin_queue.push( initial_bin );

  result = initial_bin.getResult();
  absolute_error = initial_bin.getAbsoluteError();
  int round_off_left = (T)0;
  int round_off_right = (T)0;
  
  for( int last = 1; last < d_subinterval_limit; last++ )
  {
    T result_asc_left = (T)0, result_asc_right = (T)0;
    QuadratureBin<T> problem_bin;
    QuadratureBin<T> left_half_problem_bin, right_half_problem_bin;

    // Pop bin with highest error from queue
    problem_bin = bin_queue.top();
    bin_queue.pop();

    this->bisectAndIntegrateBinInterval<Points>( 
      integrand, 
      problem_bin,
      left_half_problem_bin,
      right_half_problem_bin,
      result_asc_left,
      result_asc_right );
    
    TEST_FOR_EXCEPTION( this->subintervalTooSmall<Points>(
                                      left_half_problem_bin.getLowerLimit()
                                      right_half_problem_bin.getLowerLimit(),
                                      right_half_problem_bin.getUpperLimit() ),
                        Utility::IntegratorException,
                        "Error: Maximum number of subdivisions reached before "
                        "convergence!" );

    bin_queue.push( left_half_problem_bin );
    bin_queue.push( right_half_problem_bin );

    // Improve previous approximations to integral and error and test for
    // accuracy
    absolute_error += left_half_problem_bin.getAbsoluteError() +
      right_half_problem_bin.getAbsoluteError() -
      problem_bin.getAbsoluteError();
    
    result += left_half_problem_bin.getResult() +
      right_half_problem_bin.getResult() - problem_bin.getResult();

    // Check that the roundoff error is not too high
    this->checkRoundoffError( bin, 
                              left_half_problem_bin, 
                              right_half_problem_bin,    
                              result_asc_left,
                              result_asc_right,
                              round_off_left,
                              round_off_right,
                              last+1 );

    tolerance = 
      this->getMax( d_absolute_error_tol,
                    d_relative_error_tol*fabs( result ) );

    if( absolute_error <= tolerance )
      break;

    if( last+1 == d_subinterval_limit )
    {
      THROW_EXCEPTION( Utility::IntegratorException,
                       "Error: Maximum number of subdivisions reached before "
                       "convergence!" );
    }
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
template<typename T>
template<typename Functor>
void UnitAwareGaussKronrodIntegrator<T>::integrateAdaptivelyWynnEpsilon( 
    Functor& integrand,
    const Teuchos::ArrayView<T>& points_of_interest,
    T& result,
    T& absolute_error ) const
{
  // Check that the points of interest are in ascending order
  testPrecondition( Sort::isSortedAscending( points_of_interest.begin(),
                                             points_of_interest.end(),
                                             true ) );
  // Check that the number of points don't exceed the subinterval limit
  testPrecondition( points_of_interest.size() < d_subinterval_limit );
  // Check that there are at least two points
  testPrecondition( points_of_interest.size() > 1 );

  BinArray bin_array( d_subinterval_limit );
  Teuchos::Array<int> bin_order( number_of_intervals );

  // Initialize the bins of the integral and check for fast convergence
  {
    const bool converged;

    try{
      converged = this->initializeBinsWynnEpsilon( integrand,
                                                   points_of_interest,
                                                   bin_order,
                                                   result,
                                                   absolute_error );
    }
    EXCEPTION_CATCH_RETHROW( Utility::IntegratorException,
                             "Error: Could not initialize the intervals of "
                             "the adaptive Wynn Epsilon integration!" );

    // The integral has already converged - finish
    if( converged )
      return;
  }

  // initialize
  bin_extrapolated_result[0] = total_area;
  absolute_error = std::numeric_limits<T>::max();
  int nr_max = 0;
  int number_of_extrapolated_calls = 0;
  int number_of_extrapolated_intervals = 0;
  int ktmin = 0;
  bool extrapolate = false;
  bool no_extrapolation_allowed = false;
  bool bad_integration_behavior = false;
  T error_over_large_bins = total_error;
  T extrapolated_tolerance = tolerance;
  int max_level = 1;
  int round_off_1 = 0;
  int round_off_2 = 0;
  int round_off_3 = 0;
  int ierro = 0;
 
  T error_correction = (T)0;

  int ksgn;
  if ( fabs_total_area >= (T)1 - ((T)50*round_off) )
  {
    ksgn = 1;
  }
  else
  {
    ksgn = -1;
  }

  //number_of_intervals++;
  for ( number_of_intervals; number_of_intervals < d_subinterval_limit; number_of_intervals++ )
  {
    T area_12 = (T)0, error_12 = (T)0;
    T result_asc_1 = (T)0, result_asc_2 = (T)0;
    T smallest_bin_size = (T)0; // 1.5*smallest bin size 

    ExtrapolatedQuadratureBin<T> bin_1, bin_2;

    // Set bin to interval with largest error
    bin = bin_array[bin_order[nr_max]];

    this->bisectAndIntegrateBinInterval<21>( integrand, 
                                             bin,
                                             bin_1,
                                             bin_2,
                                             result_asc_1,
                                             result_asc_2 );
    
    bin_1.setLevel( bin.getLevel() + 1 );
    bin_2.setLevel( bin.getLevel() + 1 );

    // Improve previous approximations to integral and error and test for
    // accuracy
    total_error += bin_1.getAbsoluteError() + bin_2.getAbsoluteError() -
      bin.getAbsoluteError();
    total_area += bin_1.getResult() + bin_2.getResult() - bin.getResult();

    // Check that the roundoff error is not too high
    this->checkRoundoffError( bin, 
                              bin_1, 
                              bin_2,    
                              result_asc_1,
                              result_asc_2,
                              round_off_1,
                              round_off_2,
                              round_off_3,
                              extrapolate,
                              number_of_intervals );

    // Update and sort bin order
    this->sortBins( bin_order, 
                    bin_array,
                    bin_1,
                    bin_2,
                    number_of_intervals,
                    nr_max );

    tolerance = 
      this->getMax( d_absolute_error_tol,
                    d_relative_error_tol*fabs( total_area ) );

    if ( total_error <= tolerance )
      break;

    TEST_FOR_EXCEPTION( number_of_intervals+1 == d_subinterval_limit, 
                        Utility::IntegratorException,
                        "Maximum number of subdivisions reached" );

    TEST_FOR_EXCEPTION( subintervalTooSmall<21>( bin_1.lower_limit, 
                                                 bin_2.lower_limit, 
                                                 bin_2.upper_limit ), 
                        Utility::IntegratorException,
                        "Maximum number of subdivisions reached" );

    if ( round_off_2 >= 5 )
    {
      bad_integration_behavior = true;
    }

    if( no_extrapolation_allowed )
      continue; // go to next for loop iteration without extrapolating

    error_over_large_bins -= bin.getAbsoluteError();

    if( bin_1.getLevel() + 1 <= max_level )
    {
      error_over_large_bins += error_12;
    }

    bin = bin_array[bin_order[nr_max]];
    // Test whether the interval to be bisected next is the smallest interval
    if( !extrapolate )
    {
      if( bin.getLevel() + 1 <= max_level )
        continue; // go to next for loop iteration without extrapolating

      extrapolate = true;
      nr_max = 1; 
    }

    if( bad_integration_behavior != true && 
        error_over_large_bins > extrapolated_tolerance )
    {
      int id = nr_max;
      int size = number_of_intervals;

      if( number_of_intervals >  2 + d_subinterval_limit/2 )
        size = d_subinterval_limit + 3 - number_of_intervals;

      bool still_have_large_bins = false;
      for( int k = id; k < size; k++ )
      {
        bin = bin_array[bin_order[nr_max]];
        if( bin.getLevel() + 1 <= max_level )
        {
          still_have_large_bins = true;
          continue;
        }   
        nr_max++;    
      }
    
      if( still_have_large_bins )
        continue;
    }

    // Perform extrapolation
    T extrapolated_result = (T)0;
    T extrapolated_error = (T)0;

    number_of_extrapolated_intervals++;
    bin_extrapolated_result[number_of_extrapolated_intervals] = total_area;

    if( number_of_extrapolated_intervals < 2 )
    {
      nr_max = 0;
      extrapolate = false;
      max_level++;
      error_over_large_bins = total_error;
      continue;
    }

    Teuchos::Array<T> last_three_results( 3 );

    this->getWynnEpsilonAlgorithmExtrapolation( 
                                              bin_extrapolated_result,
                                              last_three_results,
                                              extrapolated_result, 
                                              extrapolated_error,  
                                              number_of_extrapolated_intervals,
                                              number_of_extrapolated_calls );

    ktmin++;

    TEST_FOR_EXCEPTION( ktmin > 5 &&
                        absolute_error < (1/(T)1000)*total_error, 
                        Utility::IntegratorException,
                        "Error: The integral is probably divergent "
                        "(or slowly convergent)!" );

    if( extrapolated_error < absolute_error )
    {
      ktmin = 0;
      absolute_error = extrapolated_error;
      result = extrapolated_result;
      error_correction = error_over_large_bins;
      extrapolated_tolerance = 
        this->getMax( d_absolute_error_tol, 
                      d_relative_error_tol*fabs( extrapolated_result ) ); 

      if( absolute_error <= extrapolated_tolerance )
        break;
    }

    // Prepare bisection of the smallest interval.
    if( number_of_extrapolated_intervals == 0 )
      no_extrapolation_allowed = true;

    nr_max = 0;
    extrapolate = false;
    max_level++;
    error_over_large_bins = total_error;
  } // end main for loop

  //  Set final result and error estimate.
  if( absolute_error == std::numeric_limits<T>::max() )
  {
    //  Compute global integral sum.  
    T long_result = (T)0;
    Teuchos::Array<int>::reverse_iterator j =  bin_order.rbegin();
    // Sum result over all bins
    for( j; j != bin_order.rend(); j++ )
    {
      bin = bin_array[*j];
      long_result += bin.result;
    }
    result = long_result;
    absolute_error = total_error;

    return;
  }

  if( bad_integration_behavior ) 
  {
    absolute_error += error_correction;

    TEST_FOR_EXCEPTION( bad_integration_behavior, 
                        Utility::IntegratorException,
                        "Error: Extremely bad integrand behavior occurs at "
                        "some points of the integration interval!" );

  }

  // Test on divergence.
  if( ksgn == (-1) && this->getMax( fabs(result), fabs(total_area) ) <=
      (total_area_abs/(T)100) )
  {
    return;
  }

  TEST_FOR_EXCEPTION( (1/(T)100) > result/total_area ||
                      result/total_area > (T)100 ||
                      total_error > fabs( total_area ), 
                      Utility::IntegratorException,
                      "Error: The input is invalid! Because the absolute "
                      "error tol < 0 and the relative error tol < 0, the "
                      "result and absolute_error are set to zero." );

  return;
}

// Conduct the first iteration of the adapative integration
template<typename T>
template<int Points, typename Functor>
bool UnitAwareGaussKronrodIntegrator<T>::initializeBinsWynnEpsilon(
                               Functor& integrand,
                               const Teuchos::ArrayView<T>& points_of_interest,
                               BinArray& bin_array,
                               Teuchos::Array<unsigned>& bin_order,
                               T& result,
                               T& absolute_error ) const
{
  // Initialize the bin array
  bin_array.resize( d_subinterval_limit );

  // Initialize the bin order array
  int number_of_intervals = points_of_interest.size() - 1;
  bin_order.resize( number_of_intervals );

  // Keep track of which bins need their errors rescaled
  Teuchos::Array<bool> rescale_bin_error( number_of_intervals );

  // Initialize the absolute error
  absolute_error = (T)0;
  
  T total_area = (T)0, total_area_abs = (T)0;
  T total_error = (T)0;

  // Compute the integration between the points of interest
  for( int i = 0; i < number_of_intervals; i++ )
  {
    bin_array[i].setLowerLimit( points_of_interest[i] );
    bin_array[i].setUpperLimit( points_of_interest[i+1] );

    T local_result, local_absolute_error;
    T result_abs = (T)0;
    T result_asc = (T)0;

    this->integrateWithPointRule<21>( integrand, 
                                      bin_array[i].getLowerLimit()
                                      bin_array[i].getUpperLimit(), 
                                      local_result, 
                                      local_absolute_error,
                                      result_abs, 
                                      result_asc );

    bin_array[i].setResultAndAbsoluteError( local_result,
                                            local_absolute_error );
    bin_array[i].setLevel( 0 );
    
    total_area += bin_array[i].getResult();
    absolute_error += bin_array[i].getAbsoluteError();
    total_area_abs += result_abs;

    // Check if the bin error needs to be rescaled
    if( bin_array[i].getAbsoluteError() == result_asc &&
        bin_array[i].getAbsoluteError() != (T)0 )
      rescale_bin_error[i] = true;
    else
      rescale_bin_error[i] = false;
  }
  
  // Compute error approximations for integrals between the points of interest
  for( int i = 0; i < number_of_intervals; i++ )
  {
    if( rescale_bin_error[i] ) 
      bin_array[i].setAbsoluteError( absolute_error );

    total_error += bin_array[i].getAbsoluteError();

    bin_order[i] = i;
  }

  // Sort bins from greatest to lowest error
  std::sort( bin_array.rend()-number_of_intervals, bin_array.rend() );

  // Calculate the convergence parameters
  T tolerance = this->getMax( d_absolute_error_tol,
                              d_relative_error_tol*fabs( total_area ) );

  T round_off = std::numeric_limits<T>::epsilon()*total_area_abs;

  // Check for convergence
  if( absolute_error <= tolerance )
  {
    result = total_area;
    
    return true;
  }
  // Check for an insufficient number of iterations
  else if( d_subinterval_limit == 1 )
  {
    THROW_EXCEPTION( Utility::IntegratorException,
                     "Error: A maximum of one iteration was insufficient for "
                     "convergence within the tolerance!" );
  }
  // Check for bad behavior
  else if( absolute_error <= (T)100*round_off &&
           absolute_error > tolerance )
  {
    THROW_EXCEPTION( Utility::IntegratorException,
                     "Error: Cannot reach tolerance because of roundoff error"
                     " on the first attempt!" );
  }
  // Normal behavior - not converged yet
  else
    return false;
}

// Test if subinterval is too small
template<typename T>
template<int Points>
inline bool GaussKronrodIntegrator<T>::subintervalTooSmall( 
        T& lower_limit_1, 
        T& lower_limit_2, 
        T& upper_limit_2 ) const
{
  int c = Points/10;
  T max = getMax( fabs ( lower_limit_1 ), fabs ( upper_limit_2 ) );
  T epsilon = (T)1000*c*std::numeric_limits<T>::epsilon();
  T min = (T)10000*std::numeric_limits<T>::min();

  if ( max <= ( (T)1 + epsilon ) * ( fabs( lower_limit_2 ) + min ) )
    return true;
  else
    return false;
};

// Calculate the quadrature upper and lower integrand values at an abscissa
template<typename FloatType, typename ArgUnit, typename IntegrandUnit>
template<typename Functor>
void UnitAwareGaussKronrodIntegrator<T>::calculateQuadratureIntegrandValuesAtAbscissa( 
                               Functor& integrand, 
                               ArgQuantity abscissa,
                               ArgQuantity half_length,
                               ArgQuantity midpoint,
                               IntegrandQuantity& integrand_value_lower,
                               IntegrandQuantity& integrand_value_upper ) const
{
  ArgQuantity weighted_abscissa = getRawQuantity(half_length)*abscissa;

  integrand_value_lower = integrand( midpoint - weighted_abscissa );
  integrand_value_upper = integrand( midpoint + weighted_abscissa );
}; 


// Bisect and integrate the given bin interval
template<typename T>
template<int Points, typename Functor, typename Bin>
void UnitAwareGaussKronrodIntegrator<T>::bisectAndIntegrateBinInterval( 
    Functor& integrand, 
    const Bin& bin,
    Bin& bin_1,
    Bin& bin_2,
    T& bin_1_asc,
    T& bin_2_asc ) const
{
    // Bisect the bin with the largest error estimate into bin 1 and bin 2 

    // Set bin_1
    bin_1.lower_limit = bin.lower_limit;
    bin_1.upper_limit = (1/(T)2)  * ( bin.lower_limit + bin.upper_limit );

    // Set bin 2
    bin_2.lower_limit = bin_1.upper_limit;
    bin_2.upper_limit = bin.upper_limit;

    T bin_1_abs, bin_2_abs;
    // Integrate over bin 1
    integrateWithPointRule<Points>(
      integrand, 
      bin_1.lower_limit, 
      bin_1.upper_limit, 
      bin_1.result, 
      bin_1.error, 
      bin_1_abs, 
      bin_1_asc );

    // Integrate over bin 2
    integrateWithPointRule<Points>(
      integrand, 
      bin_2.lower_limit, 
      bin_2.upper_limit, 
      bin_2.result, 
      bin_2.error, 
      bin_2_abs, 
      bin_2_asc );
};

} // end Utility namespace

#endif // end UTILITY_GAUSS_KRONROD_INTEGRATOR_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_GaussKronrodIntegrator_def.hpp
//---------------------------------------------------------------------------//
