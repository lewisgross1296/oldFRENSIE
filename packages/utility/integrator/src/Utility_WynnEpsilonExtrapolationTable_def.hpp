//---------------------------------------------------------------------------//
//!
//! \file   Utility_WynnEpsilonExtrapolationTable_def.hpp
//! \author Alex Robinson
//! \brief  WynnEpsilonExtrapolationTable template definitions
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_WYNN_EPSILON_EXTRAPOLATION_TABLE_DEF_HPP
#define UTILITY_WYNN_EPSILON_EXTRAPOLATION_TABLE_DEF_HPP

// Std Lib Includes
#include <algorithm>

// FRENSIE Includes
#include "Utility_ContractException.hpp"

namespace Utility{

// Constructor
template<typename IntegralQuantity>
WynnEpsilonExtrapolationTable<IntegralQuantity>::WynnEpsilonExtrapolationTable()
  : d_integrals(),
    d_last_three_integrals(),
    d_number_of_extrapolations( 0 )
{ /* ... */ }

// Get the size of the table (number of stored integrals)
template<typename IntegralQuantity>
size_t WynnEpsilonExtrapolationTable<IntegralQuantity>::getSize() const
{
  return d_integrals.size();
}

// Get the number of appended integrals (<= table size
template<typename IntegralQuantity>
size_t WynnEpsilonExtrapolationTable<IntegralQuantity>::getNumberOfAppendedIntegrals() const
{
  return d_number_of_appended_integrals;
}

// Get the number of extrapolations that have been conducted
template<typename IntegralQuantity>
size_t WynnEpsilonExtrapolationTable<IntegralQuantity>::getNumberOfExtrapolationsConducted() const
{
  return d_number_of_extrapolations;
}
  
// Get the stored integrals
template<typename IntegralQuantity>
const std::vector<IntegralQuantity>&
WynnEpsilonExtrapolationTable<IntegralQuantity>::getIntegrals() const
{
  return d_integrals;
}

// Append an integral to the table
template<typename IntegralQuantity>
void WynnEpsilonExtrapolationTable<IntegralQuantity>::appendIntegral(
                                             const IntegralQuantity& integral )
{
  if( d_integrals.size() > d_number_of_appended_integrals )
    d_integrals[d_number_of_appended_integrals] = integral;
  else
    d_integrals.push_back( integral );

  ++d_number_of_appended_integrals;
}

// Use the table to extrapolate the integral and absolute error
/*! \details
 */
template<typename IntegralQuantity>
void WynnEpsilonExtrapolationTable<IntegralQuantity>::extrapolate(
                                             IntegralQuantity& integral,
                                             IntegralQuantity& absolute_error )
{
  // Initialize the integral and absolute error
  integral = d_integrals.back();
  absolute_error = IQT::max();

  // Keep track of convergence
  bool converged = false;
  
  // Check if the table is large enough to conduct an extrapolation
  if( d_integrals.size() >= 2 )
  {
    // The original number of extrapolations conducted
    const size_t orig_num_extrapolations = d_number_of_extrapolations;
    
    // Index of the last appended integral
    const size_t n = d_number_of_appended_integrals-1;
    const size_t n_final = n;

    // Make sure the integrals array has the correct size
    while( d_integrals.size() < d_number_of_appended_integrals+2 )
      d_integrals.push_back( IQT::zero() );

    // Initialize the end of the integrals array
    d_integrals[n+2] = d_integrals[n];

    // Extrapolate and update the table
    converged =
      this->extrapolateAndUpdate( integral, absolute_error, n, n_final );

    // Shift the table if convergence has not occurred yet
    if( !converged )
    {
      this->shift( n, n_final );

      this->updateErrorUsingLastThreeExtrapolations( integral, absolute_error);
    }
  }  
}

// Extrapolate and update the table
template<typename IntegralQuantity>
bool WynnEpsilonExtrapolationTable<IntegralQuantity>::extrapolateAndUpdate(
                                              IntegralQuantity& integral,
                                              IntegralQuantity& absolute_error,
                                              const size_t n,
                                              size_t& n_final )
{
  // Keep track of convergence
  bool convergence = false;

  // Half of the index of the last appended integral
  const size_t half_n = n/2;
  
  for( size_t i = 0; i < half_n; ++i )
  {
    // Get the integrals of interest
    IntegralQuantity possible_integral = d_integrals[n-2*i+2];
      
    IntegralQuantity integral_0 = d_integrals[n-2*i-2];
    IntegralQuantity integral_1 = d_integrals[n-2*i-1];
    IntegralQuantity integral_2 = possible_integral;
    
    // Compute the error estimates and the tolerances
    IntegralQuantity integral_0_abs = fabs(integral_0);
    IntegralQuantity integral_1_abs = fabs(integral_1);
    IntegralQuantity integral_2_abs = fabs(integral_2);
    
    IntegralQuantity delta_21 = integral_2 - integral_1;
    IntegralQuantity delta_10 = integral_1 - integral_0;
    
    IntegralQuantity error_21 = fabs(delta_21);
    IntegralQuantity tol_21 =
      std::max( integral_2_abs, integral_1_abs )*QT::epsilon();
    IntegralQuantity error_10 = fabs(delta_10);
    IntegralQuantity tol_10 =
      std::max( integral_1_abs, integral_0_abs )*QT::epsilon();
    
    // Check if the convergence has been achieved (integral 1, 2, 3 are
    // all within machine precision).
    if( error_21 <= tol_21 && error_10 <= tol_10 )
    {
      integral = possible_integral;
      
      IntegralQuantity possible_absolute_error_a = error_21 + error_10;
      IntegralQuantity possible_absolute_error_b =
        5*QT::epsilon()*fabs(possible_integral);
      
      absolute_error = std::max( possible_absolute_error_a,
                                 possible_absolute_error_b );
      
      converged = true;
      
      break;
    }
      
    // Convergence has not been achieved - update the table
    else
    {
      IntegralQuantity integral_3 = d_integrals[n-2*i];
      d_integrals[n-2*i] = integral_1;
      
      IntegralQuantity integral_3_abs = fabs(integral_3_abs);

      IntegralQuantity delta_13 = integral_1 - integral_3;
      IntegralQuantity error_13 = fabs(delta_13);
      IntegralQuantity tol_13 =
        std::max( integral_1_abs, integral_3_abs )*QT::epsilon();

      // Omit elements that are close together (they will be overwritten)
      if( error_13 <= tol_13 || error_21 <= tol_21 || error_10 <= tol_10 )
      {
        n_final = 2*i;

        break;
      }

      // Check for irregular behavior in the table and omit the problem element
      typename IQT::RawType behavior_test_value =
        integral_1*(1/delta_13 + 1/delta_21 - 1/delta_10);

      if( behavior_test_value <= 0.0001 )
      {
        n_final = 2*i;

        break;
      }

      // Regular behavior - update the table values
      else
      {
        possible_integral = integral_1*( 1 + 1/behavior_test_value );
        
        d_integrals[n-2*i] = possible_integral;

        const IntegralQuantity possible_absolute_error =
          error_21 + fabs(possible_integral - integral_2) + error_10;

        // Only update the integral and absolute error if the extrapolation
        // improved their estimates
        if( possible_absolute_error < absolute_error )
        {
          integral = possible_integral;
          
          absolute_error = possible_absolute_error;
        }
      }
    }
  }
}

// Shift the table
template<typename IntegralQuantity>
void WynnEpsilonExtrapolationTable<IntegralQuantity>::shiftTable(
                                                         const size_t n,
                                                         const size_t n_final )
{
  // Shift back 1 element
  {
    const size_t n_limit = 49;

    if( n_final == n_limit )
      n_final = 2*(n_limit/2);
  }

  // Half of the index of the last appended integral
  const size_t half_n = n/2;

  // Shift to the left by two elements (starting with 2nd element)
  if( n%2 == 1 )
  {
    for( size_t i = 0; i < half_n; ++i )
    {
      if( d_integrals.size() > 3+i*2 )
        d_integrals[1+i*2] = d_integrals[3+i*2];
      else
        d_integrals[1+i*2] = IQT::zero();
    }

    // Shift to the left by two elements (starting with 1st element)
    else
    {
      for( size_t i = 0; i < half_n; ++i )
      {
        if( d_integrals.size() > i*2 + 2)
          d_integrals[i*2] = d_integrals[i*2+2];
        else
          d_integrals[i*2] = IQT::zero();
      }
    }
  }

  // Shift by the difference between n and n_final
  if( n != n_final )
  {
    for( size_t i = 0; i < n_final; ++i )
    {
      d_integrals[i] = d_integrals[n - n_final + i];
    }
  }

  // Set the number of appended integrals counter
  d_number_of_appended_integrals = n_final + 1;
}

// Update the error estimate using the last three extrapolated integrals
template<typename IntegralQuantity>
void WynnEpsilonExtrapolationTable<IntegralQuantity>::updateErrorUsingLastThreeExtrapolations(
                                             const IntegralQuantity& integral,
                                             IntegralQuantity& absolute_error )
{
  // Store the last three extrapolated integrals
  if( d_number_of_extrapolations < 3 )
  {
    d_last_three_integrals[d_number_of_extrapolations] = integral;

    absolute_error = IQT::max();
  }
  else
  {
    // Estimate the absolute error using the last three extrapolated integrals
    absolute_error = fabs(integral - d_last_three_integrals[2]) +
      fabs(integral - d_last_three_integrals[1]) +
      fabs(integral - d_last_three_integrals[0]);
      
    d_last_three_integrals[0] = d_last_three_integrals[1];
    d_last_three_integrals[1] = d_last_three_integrals[2];
    d_last_three_integrals[2] = integral;
  }

  // Set the absolute error
  absolute_error = std::max( absolute_error, 5*QT::epsilon()*fabs(integral) );

  // Increment the number of extrapolations conducted
  ++d_number_of_extrapolations;
}

} // end Utility namespace

#endif // end Utility_WYNN_EPSILON_EXTRAPOLATION_TABLE_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_WynnEpsilonExtrapolationTable_def.hpp
//---------------------------------------------------------------------------//
