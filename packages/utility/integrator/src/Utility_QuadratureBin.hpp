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
template<typename IndepQuantity, typename ResultQuantity = IndepQuantity>
class QuadratureBin
{

protected:

  // Quantity traits for the IndepQuantity
  typedef Utility::QuantityTraits<IndepQuantity> IQT;

  // Quantity traits for the ResultQuantity
  typedef Utility::QuantityTraits<DepQuantity> RQT;

public:

  //! Default constructor
  QuadratureBin();

  //! Constructor
  QuadratureBin( const IndepQuantity& lower_limit,
                 const IndepQuantity& upper_limit );

  //! Copy constructor
  QuadratureBin( const QuadratureBin& other_bin );

  //! Assignment operator
  QuadratureBin& operator=( const QuadratureBin& other_bin );

  //! Comparison operator
  bool operator<( const QuadratureBin& other_bin );

  //! Destructor
  virtual ~QuadratureBin()
  { /* ... */ }

  //! Get the lower limit of the bin
  const IndepQuantity& getLowerLimit() const;

  //! Get the upper limit of the bin
  const IndepQuantity& getUpperLimit() const;

  //! Set the integral of the bin
  void setIntegral( const ResultQuantity& integral );

  //! Get the integral of the bin
  const ResultQuantity& getIntegral() const;

  //! Get the integral of the bin
  ResultQuantity& getIntegral();

  //! Set the integral abs of the bin
  void setIntegralAbs( const ResultQuantity& integral_abs );

  //! Get the integral abs of the bin
  const ResultQuantity& getIntegralAbs() const;

  //! Get the integral abs of the bin
  ResultQuantity& getIntegralAbs();

  //! Set the integral asc of the bin
  void setIntegralAsc( const ResultQuantity& integral_asc );

  //! Get the integral asc of the bin
  const ResultQuantity& getIntegralAsc() const;

  //! Get the integral asc of the bin
  ResultQuantity& getIntegralAsc();

  //! Set the error of the bin integral
  void setAbsoluteError( const ResultQuantity& error );

  //! Get the error of the integral
  const ResultQuantity& getAbsoluteError() const;

  //! Get the error of the integral
  ResultQuantity& getAbsoluteError();

  //! Set the integral and error of the bin
  void setIntegralData( const ResultQuantity& integral,
                        const ResultQuantity& integral_abs,
                        const ResultQuantity& integral_asc,
                        const ResultQuantity& error );

private:

  // The lower limit of the bin
  IndepQuantity d_lower_limit;

  // The upper limit of the bin
  IndepQuantity d_upper_limit;

  // The integral of the bin
  ResultQuantity d_integral;

  // The integral abs of the bin
  ResultQuantity d_integral_abs;

  // The integral asc of the bin
  ResultQuantity d_integral_asc;

  // The absolute error of the bin
  ResultQuantity d_absolute_error;
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
