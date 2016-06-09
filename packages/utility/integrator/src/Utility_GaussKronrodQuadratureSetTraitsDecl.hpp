//---------------------------------------------------------------------------//
//!
//! \file   Utility_GaussKronrodQuadratureSetTraitsDecl.hpp
//! \author Luke Kersting
//! \brief  Gauss-Kronrod quadrature set traits declaration
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_GAUSS_KRONROD_QUADRATURE_SET_TRAITS_DECL_HPP
#define UTILITY_GAUSS_KRONROD_QUADRATURE_SET_TRAITS_DECL_HPP

// Std Lib Includes
#include <vector>

// FRENSIE Includes
#include "Utility_UndefinedTraits.hpp"
#include "Utility_UnitTraits.hpp"

/*! \defgroup gauss_kronrod_quad_traits Gauss-Kronrod Quadrature Set Traits
 * \ingroup traits
 *
 * Gauss-Kronrod quadrature sets are defined for specific point rules. If
 * an invalid point rule is requested, the default quadrature set traits
 * class will not compile. The compile time error message is defined by
 * the Utility::UndefinedTraits struct.
 */

namespace Utility{

  /*! Gauss-Kronrod quadrature set traits 
   * \ingroup gauss_kronrod_quad_traits
   */
  template<int Points, typename Unit = void, typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits
  {
    //! The abscissae quantity
    typedef typename UnitTraits<Unit>::template GetQuantityType<FloatType>::type AbscissaQuantity;

    //! The weight quantity (integrand q*weight q == integral q)
    typedef AbscissaQuantity WeightQuantity;

    //! Get the Gauss quadrature weights 
    static inline const std::vector<WeightQuantity> getGaussWeights()
    { (void)UndefinedTraits<int>::notDefined(); return std::vector<FloatType>(); }
    
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<WeightQuantity> getKronrodWeights()
    { (void)UndefinedTraits<int>::notDefined(); return std::vector<FloatType>(); }

    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<AbscissaeQuantity> getKronrodAbscissae()
    { (void)UndefinedTraits<int>::notDefined(); return std::vector<AbscissaQuantity>(); }
  };

} // end Utility namespace

#endif // end UTILITY_GAUSS_KRONROD_QUADRATURE_SET_TRAITS_DECL_HPP

//---------------------------------------------------------------------------//
// end Utility_GaussKronrodQuadratureSetTraitsDecl.hpp
//---------------------------------------------------------------------------//
