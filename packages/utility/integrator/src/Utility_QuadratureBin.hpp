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
  QuadratureBin( const IndepQuantity lower_limit,
                 const IndepQuantity upper_limit );

  //! Assignment operator
  QuadratureBin& operator=( const QuadratureBin& other_bin );

  //! Comparison operator
  bool operator<( const QuadratureBin& other_bin );

  //! Destructor
  virtual ~QuadratureBin()
  { /* ... */ }

  //! Get the lower limit of the bin
  IndepQuantity getLowerLimit() const;

  //! Get the upper limit of the bin
  IndepQuantity getUpperLimit() const;

  //! Set the result of the bin
  void setResult( const ResultQuantity result );

  //! Get the result of the bin
  ResultQuantity getResult() const;

  //! Set the error of the result
  void setAbsoluteError( const ResultQuantity error );

  //! Get the error of the result
  ResultQuantity getAbsoluteError() const;

  //! Set the result and error of the bin
  void setResultAndAbsoluteError( const ResultQuantity result,
                                  const ResultQuantity error );

private:

  // The lower limit of the bin
  IndepQuantity d_lower_limit;

  // The upper limit of the bin
  IndepQuantity d_upper_limit;

  // The result of the bin
  ResultQuantity d_result;

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
