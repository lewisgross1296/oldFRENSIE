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
  template<typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<15,FloatType>
  {
    //! Get the Gauss quadrature weights 
    static inline const std::vector<FloatType>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<FloatType>& getKronrodWeights();
      
    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<FloatType>& getKronrodAbscissae();
        
    private:

    // Initialize the gauss weight array
    static std::vector<FloatType> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<FloatType> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<FloatType> initializeKronrodAbscissae();
  };

  /*! Gauss-Kronrod quadrature set traits 21 point rule
   * \ingroup gauss_kronrod_quad_traits
   */
  template<typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<21,FloatType>
  {
    //! Get the Gauss quadrature weights 
    static inline const std::vector<FloatType>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<FloatType>& getKronrodWeights();

    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<FloatType>& getKronrodAbscissae();

    private:

    // Initialize the gauss weight array
    static std::vector<FloatType> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<FloatType> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<FloatType> initializeKronrodAbscissae();
  };

  /*! Gauss-Kronrod quadrature set traits 31 point rule
   * \ingroup gauss_kronrod_quad_traits
   */
  template<typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<31,FloatType>
  {
    //! Get the Gauss quadrature weights 
    static inline const std::vector<FloatType>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<FloatType>& getKronrodWeights();

    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<FloatType>& getKronrodAbscissae();

    private:

    // Initialize the gauss weight array
    static std::vector<FloatType> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<FloatType> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<FloatType> initializeKronrodAbscissae();
  };

  /*! Gauss-Kronrod quadrature set traits 41 point rule
   * \ingroup gauss_kronrod_quad_traits
   */
  template<typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<41,FloatType>
  {
    //! Get the Gauss quadrature weights 
    static inline const std::vector<FloatType>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<FloatType>& getKronrodWeights();

    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<FloatType>& getKronrodAbscissae();
    
    private:

    // Initialize the gauss weight array
    static std::vector<FloatType> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<FloatType> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<FloatType> initializeKronrodAbscissae();
  };

  /*! Gauss-Kronrod quadrature set traits 51 point rule
   * \ingroup gauss_kronrod_quad_traits
   */
  template<typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<51,FloatType>
  {
    //! Get the Gauss quadrature weights 
    static inline const std::vector<FloatType>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<FloatType>& getKronrodWeights();

    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<FloatType>& getKronrodAbscissae();

    private:

    // Initialize the gauss weight array
    static std::vector<FloatType> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<FloatType> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<FloatType> initializeKronrodAbscissae();
  };

  /*! Gauss-Kronrod quadrature set traits 61 point rule
   * \ingroup gauss_kronrod_quad_traits
   */
  template<typename FloatType = double>
  struct GaussKronrodQuadratureSetTraits<61,FloatType>
  {
    //! Get the Gauss quadrature weights 
    static inline const std::vector<FloatType>& getGaussWeights();
        
    //! Get the Kronrod quadrature weights 
    static inline const std::vector<FloatType>& getKronrodWeights();

    //! Get the Kronrod quadrature abscissae
    static inline const std::vector<FloatType>& getKronrodAbscissae();

    private:

    // Initialize the gauss weight array
    static std::vector<FloatType> initializeGaussWeights();

    // Initialize the kronrod weight array
    static std::vector<FloatType> initializeKronrodWeights();

    // Initialize the kronrod abscissae array
    static std::vector<FloatType> initializeKronrodAbscissae();
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
