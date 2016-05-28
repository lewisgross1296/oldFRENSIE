//---------------------------------------------------------------------------//
//!
//! \file   Utility_ExtrapolatedQuadratureBin.hpp
//! \author Alex Robinson
//! \brief  The extrapolated quadrature bin definition
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_EXTRAPOLATED_QUADRATURE_BIN_DEF_HPP
#define UTILITY_EXTRAPOLATED_QUADRATURE_BIN_DEF_HPP

// FRENSIE Includes
#include "Utility_ContractException.hpp"

namespace Utility{

//! Default constructor
template<typename IndepQuantity, typename ResultQuantity = IndepQuantity>
ExtrapolatedQuadratureBin<IndepQuantity,ResultQuantity>::ExtrapolatedQuadratureBin()
  : QuadratureBin<IndepQuantity,ResultQuantity>(),
    d_level( 0u )
{ /* ... */ }

//! Constructor
template<typename IndepQuantity, typename ResultQuantity = IndepQuantity>
ExtrapolatedQuadratureBin<IndepQuantity,ResultQuantity>::ExtrapolatedQuadratureBin(
                                              const IndepQuantity lower_limit,
                                              const IndepQuantity upper_limit )
  : QuadratureBin<IndepQuantity,ResultQuantity>( lower_limit, upper_limit ),
    d_level( 0u )
{ /* ... */ }

//! Assignment operator
template<typename IndepQuantity, typename ResultQuantity = IndepQuantity>
ExtrapolatedQuadratureBin<IndepQuantity,ResultQuantity>&
ExtrapolatedQuadratureBin<IndepQuantity,ResultQuantity>::operator=(
     const ExtrapolatedQuadratureBin<IndepQuantity,ResultQuantity>& other_bin )
{
  QuadratureBin<IndepQuantity,ResultQuantity>::operator=( other_bin );

  if( this != &other_bin )
    d_level = other_bin.d_level;
}

//! Set the bin level
template<typename IndepQuantity, typename ResultQuantity = IndepQuantity>
void ExtrapolatedQuadratureBin<IndepQuantity,ResultQuantity>::setLevel(
                                                         const unsigned level )
{
  d_level = level;
}

//! Get the bin level
template<typename IndepQuantity, typename ResultQuantity = IndepQuantity>
unsigned
ExtrapolatedQuadratureBin<IndepQuantity,ResultQuantity>::getLevel() const
{
  return d_level;
}
  
} // end Utility namespace

#endif // end UTILITY_EXTRAPOLATED_QUADRATURE_BIN_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_ExtrapolatedQuadratureBin_def.hpp
//---------------------------------------------------------------------------//
