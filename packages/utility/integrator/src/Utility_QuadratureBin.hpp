//---------------------------------------------------------------------------//
//!
//! \file   Utility_QuadratureBin.hpp
//! \author Alex Robinson, Luke Kersting
//! \brief  The quadrature bin declaration (previously BinTraits)
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_QUADRATURE_BIN_HPP
#define UTILITY_QUADRATURE_BIN_HPP

// FRENSIE Includes
#include "Utility_QuantityTraits.hpp"

namespace Utility{

//! The quadrature bin class
template<typename ArgQuantity, typename IntegralQuantity = ArgQuantity>
class QuadratureBin
{

protected:

  //! Quantity traits for the ArgQuantity
  typedef Utility::QuantityTraits<ArgQuantity> AQT;

  //! Quantity traits for the IntegralQuantity
  typedef Utility::QuantityTraits<IntegralQuantity> IQT;

public:

  //! Default constructor
  QuadratureBin();

  //! Constructor
  QuadratureBin( const ArgQuantity& lower_limit,
                 const ArgQuantity& upper_limit );

  //! Copy constructor
  QuadratureBin( const QuadratureBin& other_bin );

  //! Assignment operator
  QuadratureBin& operator=( const QuadratureBin& other_bin );

  //! Destructor
  virtual ~QuadratureBin()
  { /* ... */ }

  //! Get the lower limit of the bin
  const ArgQuantity& getLowerLimit() const;

  //! Get the upper limit of the bin
  const ArgQuantity& getUpperLimit() const;

  //! Set the integral of the bin
  void setIntegral( const IntegralQuantity& integral );

  //! Get the integral of the bin
  const IntegralQuantity& getIntegral() const;

  //! Get the integral of the bin
  IntegralQuantity& getIntegral();

  //! Set the integral abs of the bin
  void setIntegralAbs( const IntegralQuantity& integral_abs );

  //! Get the integral abs of the bin
  const IntegralQuantity& getIntegralAbs() const;

  //! Get the integral abs of the bin
  IntegralQuantity& getIntegralAbs();

  //! Set the integral asc of the bin
  void setIntegralAsc( const IntegralQuantity& integral_asc );

  //! Get the integral asc of the bin
  const IntegralQuantity& getIntegralAsc() const;

  //! Get the integral asc of the bin
  IntegralQuantity& getIntegralAsc();

  //! Set the error of the bin integral
  void setAbsoluteError( const IntegralQuantity& error );

  //! Get the error of the integral
  const IntegralQuantity& getAbsoluteError() const;

  //! Get the error of the integral
  IntegralQuantity& getAbsoluteError();

  //! Set the integral and error of the bin
  void setIntegralData( const IntegralQuantity& integral,
                        const IntegralQuantity& integral_abs,
                        const IntegralQuantity& integral_asc,
                        const IntegralQuantity& error );

  //! Bisect the bin
  void bisect( QuadratureBin& left_half_bin,
               QuadratureBin& right_half_bin ) const;

  //! Get the bisection level of the bin
  unsigned getBisectionLevel() const;

private:

  // Set the bin limits
  void setLimits( const ArgQuantity& lower_limit,
                  const ArgQuantity& upper_limit );

  // Increment the bisection level
  void incrementBisectionLevel();

  // Reset the integral and error of the bin
  void resetIntegralData();

  // The lower limit of the bin
  ArgQuantity d_lower_limit;

  // The upper limit of the bin
  ArgQuantity d_upper_limit;

  // The integral of the bin
  IntegralQuantity d_integral;

  // The integral abs of the bin
  IntegralQuantity d_integral_abs;

  // The integral asc of the bin
  IntegralQuantity d_integral_asc;

  // The absolute error of the bin
  IntegralQuantity d_absolute_error;

  // The bisection level (number of bisections that have been done to get
  // to this bins limits)
  unsigned d_bisection_level;
};
  
} // end Utility namespace

//---------------------------------------------------------------------------//
// Template Includes
//---------------------------------------------------------------------------//

#include "Utility_GaussKronrodQuadratureBin_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_QUADRATURE_BIN_HPP

//---------------------------------------------------------------------------//
// end Utility_QuadratureBin.hpp
//---------------------------------------------------------------------------//
