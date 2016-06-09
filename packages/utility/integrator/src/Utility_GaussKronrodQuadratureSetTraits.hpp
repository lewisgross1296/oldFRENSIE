//---------------------------------------------------------------------------//
//!
//! \file   Utility_GaussKronrodQuadratureSetTraits.hpp
//! \author Luke Kersting
//! \brief  Gauss-Kronrod quadrature set traits
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_GAUSS_KRONROD_QUADRATURE_SET_TRAITS_HPP
#define UTILITY_GAUSS_KRONROD_QUADRATURE_SET_TRAITS_HPP

// FRENSIE Includes
#include "Utility_GaussKronrodQuadratureSetTraitsDecl.hpp"

namespace Utility{

  /*! Gauss-Kronrod quadrature set traits 15 point rule
   * \ingroup gauss_kronrod_quad_traits
   */
  template<typename Unit = void, typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<15,Unit,FloatType>
  {
    //! The abscissae quantity
    typedef typename UnitTraits<Unit>::template GetQuantityType<FloatType>::type AbscissaQuantity;

    //! The weight quantity (integrand q*weight q == integral q)
    typedef AbscissaQuantity WeightQuantity;
    
    //! Get the Gauss quadrature weights 
    static inline const std::vector<WeightQuantity>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<WeightQuantity>& getKronrodWeights();
      
    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<AbscissaQuantity>& getKronrodAbscissae();
        
    private:

    // Initialize the gauss weight array
    static std::vector<WeightQuantity> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<WeightQuantity> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<AbscissaQuantity> initializeKronrodAbscissae();
  };

  /*! Gauss-Kronrod quadrature set traits 21 point rule
   * \ingroup gauss_kronrod_quad_traits
   */
  template<typename Unit = void, typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<21,Unit,FloatType>
  {
    //! The abscissa quantity
    typedef typename UnitTraits<Unit>::template GetQuantityType<FloatType>::type AbscissaQuantity;

    //! The weight quantity (integrand q*weight q == integral q)
    typedef AbscissaQuantity WeightQuantity;
    
    //! Get the Gauss quadrature weights 
    static inline const std::vector<WeightQuantity>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<WeightQuantity>& getKronrodWeights();

    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<AbscissaQuantity>& getKronrodAbscissae();

    private:

    // Initialize the gauss weight array
    static std::vector<WeightQuantity> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<WeightQuantity> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<AbscissaQuantity> initializeKronrodAbscissae();
  };

  /*! Gauss-Kronrod quadrature set traits 31 point rule
   * \ingroup gauss_kronrod_quad_traits
   */
  template<typename Unit = void, typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<31,Unit,FloatType>
  {
    //! The abscissae quantity
    typedef typename UnitTraits<Unit>::template GetQuantityType<FloatType>::type AbscissaQuantity;

    //! The weight quantity (integrand q*weight q == integral q)
    typedef AbscissaQuantity WeightQuantity;
    
    //! Get the Gauss quadrature weights 
    static inline const std::vector<WeightQuantity>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<WeightQuantity>& getKronrodWeights();

    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<AbscissaQuantity>& getKronrodAbscissae();

    private:

    // Initialize the gauss weight array
    static std::vector<WeightQuantity> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<WeightQuantity> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<AbscissaQuantity> initializeKronrodAbscissae();
  };

  /*! Gauss-Kronrod quadrature set traits 41 point rule
   * \ingroup gauss_kronrod_quad_traits
   */
  template<typename Unit = void, typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<41,Unit,FloatType>
  {
    //! The abscissae quantity
    typedef typename UnitTraits<Unit>::template GetQuantityType<FloatType>::type AbscissaQuantity;

    //! The weight quantity (integrand q*weight q == integral q)
    typedef AbscissaQuantity WeightQuantity;
    
    //! Get the Gauss quadrature weights 
    static inline const std::vector<WeightQuantity>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<WeightQuantity>& getKronrodWeights();

    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<AbscissaQuantity>& getKronrodAbscissae();
    
    private:

    // Initialize the gauss weight array
    static std::vector<WeightQuantity> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<WeightQuantity> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<AbscissaQuantity> initializeKronrodAbscissae();
  };

  /*! Gauss-Kronrod quadrature set traits 51 point rule
   * \ingroup gauss_kronrod_quad_traits
   */
  template<typename Unit = void, typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<51,Unit,FloatType>
  {
    //! The abscissae quantity
    typedef typename UnitTraits<Unit>::template GetQuantityType<FloatType>::type AbscissaQuantity;

    //! The weight quantity (integrand q*weight q == integral q)
    typedef AbscissaQuantity WeightQuantity;
    
    //! Get the Gauss quadrature weights 
    static inline const std::vector<WeightQuantity>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<WeightQuantity>& getKronrodWeights();

    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<AbscissaQuantity>& getKronrodAbscissae();

    private:

    // Initialize the gauss weight array
    static std::vector<WeightQuantity> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<WeightQuantity> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<AbscissaQuantity> initializeKronrodAbscissae();
  };

  /*! Gauss-Kronrod quadrature set traits 61 point rule
   * \ingroup gauss_kronrod_quad_traits
   */
  template<typename Unit = void, typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<61,Unit,FloatType>
  {
    //! The abscissae quantity
    typedef typename UnitTraits<Unit>::template GetQuantityType<FloatType>::type AbscissaQuantity;

    //! The weight quantity (integrand q*weight q == integral q)
    typedef AbscissaQuantity WeightQuantity;
    
    //! Get the Gauss quadrature weights 
    static inline const std::vector<WeightQuantity>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<WeightQuantity>& getKronrodWeights();

    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<AbscissaQuantity>& getKronrodAbscissae();

    private:

    // Initialize the gauss weight array
    static std::vector<WeightQuantity> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<WeightQuantity> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<AbscissaQuantity> initializeKronrodAbscissae();
  };

} // end Utility namespace

//---------------------------------------------------------------------------//
// Template Includes
//---------------------------------------------------------------------------//

#include "Utility_GaussKronrodQuadratureSetTraits_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_GAUSS_KRONROD_QUADRATURE_SET_TRAITS_HPP

//---------------------------------------------------------------------------//
// end Utility_GaussKronrodQuadratureSetTraits.hpp
//---------------------------------------------------------------------------//
