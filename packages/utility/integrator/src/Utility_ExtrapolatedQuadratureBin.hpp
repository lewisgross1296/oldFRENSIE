//---------------------------------------------------------------------------//
//!
//! \file   Utility_ExtrapolatedQuadratureBin.hpp
//! \author Alex Robinson, Luke Kersting
//! \brief  The extrapolated quadrature bin declaration (previously
//!         ExtrpolatedBinTraits)
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_EXTRAPOLATED_QUADRATURE_BIN_HPP
#define UTILITY_EXTRAPOLATED_QUADRATURE_BIN_HPP

// FRENSIE Includes
#include "Utility_QuadratureBin.hpp"

namespace Utility{

//! The extrapolated quadrature bin class
template<typename IndepQuantity, typename ResultQuantity = IndepQuantity>
class ExtrapolatedQuadratureBin : public QuadratureBin<IndepQuantity,ResultQuantity>
{

public:

  //! Default constructor
  ExtrapolatedQuadratureBin();

  //! Constructor
  ExtrapolatedQuadratureBin( const IndepQuantity lower_limit,
                             const IndepQuantity upper_limit );

  //! Destructor
  ~ExtrapolatedQuadratureBin()
  { /* ... */ }

  //! Assignment operator
  ExtrapolatedQuadratureBin& operator=(
                                  const ExtrapolatedQuadratureBin& other_bin );

  //! Set the bin level
  void setLevel( const unsigned level );

  //! Get the bin level
  unsigned getLevel() const;

private:

  // The bin level
  unsigned d_level;
};
  
} // end Utility namespace

//---------------------------------------------------------------------------//
// Template Includes
//---------------------------------------------------------------------------//

#include "Utility_ExtrapolatedQuadratureBin_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_EXTRAPOLATED_QUADRATURE_BIN_HPP

//---------------------------------------------------------------------------//
// end Utility_ExtrapolatedQuadratureBin.hpp
//---------------------------------------------------------------------------//
