//---------------------------------------------------------------------------//
//!
//! \file   Utility_WynnEpsilonExtrapolationTable.hpp
//! \author Alex Robinson
//! \brief  WynnEpsilonExtrapolationTable declaration
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_WYNN_EPSILON_EXTRAPOLATION_TABLE_HPP
#define UTILITY_WYNN_EPSILON_EXTRAPOLATION_TABLE_HPP

// Std Lib Includes
#include <vector>

// FRENSIE Includes
#include "Utility_QuantityTraits.hpp"

namespace Utility{

//! The Wynn-Epsilon extrapolation table
template<typename IntegralQuantity>
class WynnEpsilonExtrapolationTable
{

private:

  // Quantity traits for the IntegralQuantity
  typedef Utility::QuantityTraits<IntegralQuantity> IQT;

  // Quantity traits for the raw type
  typedef Utility::QuantityTraits<typename IQT::RawType> QT;
  
public:

  //! Constructor
  WynnEpsilonExtrapolationTable();

  //! Destructor
  ~WynnEpsilonExtrapolationTable()
  { /* ... */ }

  //! Get the size of the table (number of stored integrals)
  size_t getSize() const;

  //! Get the number of appended integrals (<= table size)
  size_t getNumberOfAppendedIntegrals() const;

  //! Get the number of extrapolations that have been conducted
  size_t getNumberOfExtrapolationsConducted() const;
  
  //! Get the stored integrals
  const std::vector<IntegralQuantity>& getIntegrals() const;

  //! Append an integral to the table
  void appendIntegral( const IntegralQuantity& integral );

  //! Use the table to extrapolate the integral and absolute error
  void extrapolate( IntegralQuantity& integral,
                    IntegralQuantity& absolute_error );

private:

  // Extrapolate and update the table
  bool extrapolateAndUpdate( IntegralQuantity& integral,
                             IntegralQuantity& absolute_error,
                             const size_t n,
                             size_t& n_final );

  // Shift the table
  void shiftTable( const size_t n, const size_t n_final );

  // Update the error estimate using the last three extrapolated integrals
  void updateErrorUsingLastThreeExtrapolations(
                                            const IntegralQuantity& integral,
                                            IntegralQuantity& absolute_error );

  // The stored integrals
  std::vector<IntegralQuantity> d_integrals;

  // The last three extrapolated integrals (used for updating the error)
  IntegralQuantity d_last_three_integrals[3];

  // The number of appended integrals
  size_t d_number_of_appended_integrals;

  // The number of extrapolations conducted
  size_t d_number_of_extrapolations;
};
  
} // end Utility namespace

//---------------------------------------------------------------------------//
// Template Includes.
//---------------------------------------------------------------------------//

#include "Utility_WynnEpsilonExtrapolationTable_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_WYNN_EPSILON_EXTRAPOLATION_TABLE_HPP

//---------------------------------------------------------------------------//
// end Utility_WynnEpsilonExtrapolationTable.hpp
//---------------------------------------------------------------------------//
