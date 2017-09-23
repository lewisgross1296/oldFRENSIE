//---------------------------------------------------------------------------//
//!
//! \file   Utility_CosineInterpolationPolicy.hpp
//! \author Luke Kersting
//! \brief  Cosine interpolation policy struct declarations
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_COSINE_INTERPOLATION_POLICY_HPP
#define UTILITY_COSINE_INTERPOLATION_POLICY_HPP

// FRENSIE Includes
#include "Utility_InterpolationType.hpp"
#include "Utility_InterpolationPolicy.hpp"
#include "Utility_QuantityTraits.hpp"

namespace Utility{


// /*! \brief Helper policy struct for use with the cosine interpolation structs
//  * \ingroup policies
//  */
// template<typename ParentInterpolationType>
// class InterpolationHelper
// {
  
// private:

//   // Typedef for this type
//   typedef InterpolationHelper<ParentInterpolationType> ThisType;

// public:
  
  // //! Interpolate between two processed points
  // template<typename T>
  // static T interpolate( const T processed_indep_var_0,
	// 		const T processed_indep_var,
	// 		const T processed_dep_var_0,
	// 		const T processed_slope );

  // //! Interpolate between two processed points and return the processed value
  // template<typename T>
  // static T interpolateAndProcess( const T processed_indep_var_0,
	// 			  const T processed_indep_var,
	// 			  const T processed_dep_var_0,
	// 			  const T processed_slope );

  // //! Calculate the unit base grid length (L)
  // template<typename IndepType>
  // static typename QuantityTraits<IndepType>::RawType
  // calculateUnitBaseGridLength( const IndepType grid_lower_indep_value,
  //                              const IndepType grid_upper_indep_value );

  // //! Calculate the unit base grid length from a processed grid (L)
  // template<typename T>
  // static T calculateUnitBaseGridLengthProcessed(
  //                                   const T processed_grid_lower_indep_value,
  //                                   const T processed_grid_upper_indep_value );

  // //! Calculate the unit base independent variable (eta)
  // template<typename IndepType>
  // static typename QuantityTraits<IndepType>::RawType
  // calculateUnitBaseIndepVar(
  //        const IndepType indep_var,
  //        const IndepType indep_var_min,
  //        const typename QuantityTraits<IndepType>::RawType indep_grid_length,
  //        const double tol = 1e-3 );

  // //! Calculate the unit base independent variable (eta)
  // template<typename T>
  // static T calculateUnitBaseIndepVarProcessed( const T processed_indep_var,
  //                                              const T processed_indep_var_min,
  //                                              const T indep_grid_length,
  //                                              const double tol = 1e-3 );

  // //! Calculate the independent variable (from eta)
  // template<typename IndepType>
  // static IndepType calculateIndepVar(
  //        const typename QuantityTraits<IndepType>::RawType eta,
  //        const IndepType indep_var_min,
  //        const typename QuantityTraits<IndepType>::RawType indep_grid_length,
  //        const double tol = 1e-3 );

  // //! Calculate the processed independent variable (from eta)
  // template<typename T>
  // static T calculateProcessedIndepVar( const T eta,
  //                                      const T processed_indep_var_min,
  //                                      const T indep_grid_length );

  // //! Calculate the "fuzzy" lower bound (lower bound with roundoff tolerance)
  // template<typename T>
  // static T calculateFuzzyLowerBound( const T value, const double tol = 1e-3 );

  // //! Calculate the "fuzzy" upper bound (upper bound with roundoff tolerance)
  // template<typename T>
  // static T calculateFuzzyUpperBound( const T value, const double tol = 1e-3 );

//   //! Convert the cosine variable
//   template<typename T>
//   static T convertCosineVar( const T cosine_var );
// };

/*! \brief Policy struct for interpolating data tables that require logcos-log 
 * interpolation between evaluated points.
 * \details The dependent variable is always assumed to be an angle cosine.
 * Since the angle cosine goes below zero a direct log interpolation
 * cannot be performed on it. Instead a log interpolation will be performed
 * on the change in the angle cosine (eg: 1 - mu) instead of the cosine ( mu ).
 * The interpolated value will always be cast into the cosine before retuned.
 * \ingroup policies
 */
struct LogCosLog : public InterpolationHelper<LogCosLog>
{
  //! Independent variable processing tag
  typedef LogIndepVarProcessingTag IndepVarProcessingTag;

  //! Dependent variable processing tag
  typedef LogDepVarProcessingTag DepVarProcessingTag;

  //! Get the interpolation type
  static InterpolationType getInterpolationType();

  //! Force base class template methods to be visible
  using InterpolationHelper<LogCosLog>::interpolate;
  using InterpolationHelper<LogCosLog>::interpolateAndProcess;
  using InterpolationHelper<LogCosLog>::calculateUnitBaseIndepVar;
  using InterpolationHelper<LogCosLog>::calculateUnitBaseIndepVarProcessed;
  using InterpolationHelper<LogCosLog>::calculateIndepVar;
  using InterpolationHelper<LogCosLog>::calculateProcessedIndepVar;
  using InterpolationHelper<LogCosLog>::convertCosineVar;

  //! Interpolate between two points
  template<typename IndepType, typename DepType>
  static DepType interpolate( const IndepType indep_var_0,
                              const IndepType indep_var_1,
                              const IndepType indep_var,
                              const DepType dep_var_0,
                              const DepType dep_var_1 );

  //! Interpolate between two points using the indep variable ratio (beta)
  template<typename T, typename DepType>
  static DepType interpolate( const T beta,
                              const DepType dep_var_0,
                              const DepType dep_var_1 );

  //! Interpolate between two points and return the processed value
  template<typename IndepType, typename DepType>
  static typename QuantityTraits<DepType>::RawType
  interpolateAndProcess( const IndepType indep_var_0,
                         const IndepType indep_var_1,
                         const IndepType indep_var,
                         const DepType dep_var_0,
                         const DepType dep_var_1 );

  //! Process the independent value
  template<typename T>
  static typename QuantityTraits<T>::RawType
  processIndepVar( const T cosine_var );

  //! Process the dependent value
  template<typename T>
  static typename QuantityTraits<T>::RawType
  processDepVar( const T dep_var );

  //! Recover the processed independent value
  template<typename T>
  static T recoverProcessedIndepVar( const T processed_indep_var );

  //! Recover the processed dependent value
  template<typename T>
  static T recoverProcessedDepVar( const T processed_dep_var );

  //! Test if the independent value is in a valid range (doesn't check nan/inf)
  template<typename T>
  static bool isIndepVarInValidRange( const T indep_var );

  //! Test if the dependent value is in a valid range (doesn't check nan/inf)
  template<typename T>
  static bool isDepVarInValidRange( const T dep_var );

  //! The name of the policy
  static const std::string name();
};

/*! \brief Policy struct for interpolating data tables that require log-log
 * cosine interpolation between evaluated points.
 * \details The independent variable is always assumed to be an angle cosine.
 * Since the angle cosine goes below zero a direct log interpolation
 * cannot be performed on it. Instead a log interpolation will be performed
 * on the change in the angle cosine (eg: 1 - mu) instead of the cosine ( mu ).
 * The interpolated value will always be cast into the cosine before retuned.
 * \ingroup policies
 */
struct LogLogCos : public InterpolationHelper<LogLogCos>
{
  //! Independent variable processing tag
  typedef LogIndepVarProcessingTag IndepVarProcessingTag;

  //! Dependent variable processing tag
  typedef LogDepVarProcessingTag DepVarProcessingTag;

  //! Get the interpolation type
  static InterpolationType getInterpolationType();

  //! Force base class template methods to be visible
  using InterpolationHelper<LogLogCos>::interpolate;
  using InterpolationHelper<LogLogCos>::interpolateAndProcess;
  using InterpolationHelper<LogLogCos>::calculateUnitBaseIndepVar;
  using InterpolationHelper<LogLogCos>::calculateUnitBaseIndepVarProcessed;
  using InterpolationHelper<LogLogCos>::calculateIndepVar;
  using InterpolationHelper<LogLogCos>::calculateProcessedIndepVar;
  using InterpolationHelper<LogLogCos>::convertCosineVar;

  //! Interpolate between two points
  template<typename IndepType, typename DepType>
  static DepType interpolate( const IndepType indep_var_0,
                              const IndepType indep_var_1,
                              const IndepType indep_var,
                              const DepType dep_var_0,
                              const DepType dep_var_1 );

  //! Interpolate between two points using the indep variable ratio (beta)
  template<typename T, typename DepType>
  static DepType interpolate( const T beta,
                              const DepType dep_var_0,
                              const DepType dep_var_1 );

  //! Interpolate between two points and return the processed value
  template<typename IndepType, typename DepType>
  static typename QuantityTraits<DepType>::RawType
  interpolateAndProcess( const IndepType indep_var_0,
                         const IndepType indep_var_1,
                         const IndepType indep_var,
                         const DepType dep_var_0,
                         const DepType dep_var_1 );

  //! Process the independent value
  template<typename T>
  static typename QuantityTraits<T>::RawType
  processIndepVar( const T cosine_var );

  //! Process the dependent value
  template<typename T>
  static typename QuantityTraits<T>::RawType
  processDepVar( const T dep_var );

  //! Recover the processed independent value
  template<typename T>
  static T recoverProcessedIndepVar( const T processed_indep_var );

  //! Recover the processed dependent value
  template<typename T>
  static T recoverProcessedDepVar( const T processed_dep_var );

  //! Test if the independent value is in a valid range (doesn't check nan/inf)
  template<typename T>
  static bool isIndepVarInValidRange( const T indep_var );

  //! Test if the dependent value is in a valid range (doesn't check nan/inf)
  template<typename T>
  static bool isDepVarInValidRange( const T dep_var );

  //! The name of the policy
  static const std::string name();
};

/*! \brief Policy struct for interpolating data tables that require logcos-lin
 * interpolation between evaluated points.
 * \details The dependent variable is always assumed to be an angle cosine.
 * Since the angle cosine goes below zero a direct log interpolation
 * cannot be performed on it. Instead a log interpolation will be performed
 * on the change in the angle cosine (eg: 1 - mu) instead of the cosine ( mu ).
 * The interpolated value will always be cast into the cosine before retuned.
 * \ingroup policies
 */
struct LogCosLin : public InterpolationHelper<LogCosLin>
{
  //! Independent variable processing tag
  typedef LinIndepVarProcessingTag IndepVarProcessingTag;

  //! Dependent variable processing tag
  typedef LogDepVarProcessingTag DepVarProcessingTag;

  //! Get the interpolation type
  static InterpolationType getInterpolationType();

  //! Force base class template methods to be visible
  using InterpolationHelper<LogCosLin>::interpolate;
  using InterpolationHelper<LogCosLin>::interpolateAndProcess;
  using InterpolationHelper<LogCosLin>::calculateUnitBaseIndepVar;
  using InterpolationHelper<LogCosLin>::calculateUnitBaseIndepVarProcessed;
  using InterpolationHelper<LogCosLin>::calculateIndepVar;
  using InterpolationHelper<LogCosLin>::calculateProcessedIndepVar;
  using InterpolationHelper<LogCosLin>::convertCosineVar;

  //! Interpolate between two points
  template<typename IndepType, typename DepType>
  static DepType interpolate( const IndepType indep_var_0,
                              const IndepType indep_var_1,
                              const IndepType indep_var,
                              const DepType dep_var_0,
                              const DepType dep_var_1 );

  //! Interpolate between two points using the indep variable ratio (beta)
  template<typename T, typename DepType>
  static DepType interpolate( const T beta,
                              const DepType dep_var_0,
                              const DepType dep_var_1 );

  //! Interpolate between two points and return the processed value
  template<typename IndepType, typename DepType>
  static typename QuantityTraits<DepType>::RawType
  interpolateAndProcess( const IndepType indep_var_0,
                         const IndepType indep_var_1,
                         const IndepType indep_var,
                         const DepType dep_var_0,
                         const DepType dep_var_1 );

  //! Process the independent value
  template<typename T>
  static typename QuantityTraits<T>::RawType
  processIndepVar( const T cosine_var );

  //! Process the dependent value
  template<typename T>
  static typename QuantityTraits<T>::RawType
  processDepVar( const T dep_var );

  //! Recover the processed independent value
  template<typename T>
  static T recoverProcessedIndepVar( const T processed_indep_var );

  //! Recover the processed dependent value
  template<typename T>
  static T recoverProcessedDepVar( const T processed_dep_var );

  //! Test if the independent value is in a valid range (doesn't check nan/inf)
  template<typename T>
  static bool isIndepVarInValidRange( const T indep_var );

  //! Test if the dependent value is in a valid range (doesn't check nan/inf)
  template<typename T>
  static bool isDepVarInValidRange( const T dep_var );

  //! The name of the policy
  static const std::string name();
};

/*! \brief Policy struct for interpolating data tables that require lin-logcos
 * interpolation between evaluated points.
 * \details The independent variable is always assumed to be an angle cosine.
 * Since the angle cosine goes below zero a direct log interpolation
 * cannot be performed on it. Instead a log interpolation will be performed
 * on the change in the angle cosine (eg: 1 - mu) instead of the cosine ( mu ).
 * The interpolated value will always be cast into the cosine before retuned.
 * \ingroup policies
 */
struct LinLogCos : public InterpolationHelper<LinLogCos>
{
  //! Independent variable processing tag
  typedef LogIndepVarProcessingTag IndepVarProcessingTag;

  //! Dependent variable processing tag
  typedef LinDepVarProcessingTag DepVarProcessingTag;

  //! Get the interpolation type
  static InterpolationType getInterpolationType();

  //! Force base class template methods to be visible
  using InterpolationHelper<LinLogCos>::interpolate;
  using InterpolationHelper<LinLogCos>::interpolateAndProcess;
  using InterpolationHelper<LinLogCos>::calculateUnitBaseIndepVar;
  using InterpolationHelper<LinLogCos>::calculateUnitBaseIndepVarProcessed;
  using InterpolationHelper<LinLogCos>::calculateIndepVar;
  using InterpolationHelper<LinLogCos>::calculateProcessedIndepVar;
  using InterpolationHelper<LinLogCos>::convertCosineVar;

  //! Interpolate between two points
  template<typename IndepType, typename DepType>
  static DepType interpolate( const IndepType indep_var_0,
                              const IndepType indep_var_1,
                              const IndepType indep_var,
                              const DepType dep_var_0,
                              const DepType dep_var_1 );

  //! Interpolate between two points using the indep variable ratio (beta)
  template<typename T, typename DepType>
  static DepType interpolate( const T beta,
                              const DepType dep_var_0,
                              const DepType dep_var_1 );

  //! Interpolate between two points and return the processed value
  template<typename IndepType, typename DepType>
  static typename QuantityTraits<DepType>::RawType
  interpolateAndProcess( const IndepType indep_var_0,
                         const IndepType indep_var_1,
                         const IndepType indep_var,
                         const DepType dep_var_0,
                         const DepType dep_var_1 );

  //! Process the independent value
  template<typename T>
  static typename QuantityTraits<T>::RawType
  processIndepVar( const T cosine_var );

  //! Process the dependent value
  template<typename T>
  static typename QuantityTraits<T>::RawType
  processDepVar( const T dep_var );

  //! Recover the processed independent value
  template<typename T>
  static T recoverProcessedIndepVar( const T processed_indep_var );

  //! Recover the processed dependent value
  template<typename T>
  static T recoverProcessedDepVar( const T processed_dep_var );

  //! Test if the independent value is in a valid range (doesn't check nan/inf)
  template<typename T>
  static bool isIndepVarInValidRange( const T indep_var );

  //! Test if the dependent value is in a valid range (doesn't check nan/inf)
  template<typename T>
  static bool isDepVarInValidRange( const T dep_var );

  //! The name of the policy
  static const std::string name();
};

//! Helper class used to invert the cosine interpolation policy (LogCos-Log)
template<>
struct InverseInterpPolicy<LogCosLog>
{
  //! The inverse cosine interpolation policy
  typedef LogLogCos InterpPolicy;
};

//! Helper class used to invert the cosine interpolation policy (Log-LogCos)
template<>
struct InverseInterpPolicy<LogLogCos>
{
  //! The inverse cosine interpolation policy
  typedef LogCosLog InterpPolicy;
};

//! Helper class used to invert the cosine interpolation policy (LogCos-Lin)
template<>
struct InverseInterpPolicy<LogCosLin>
{
  //! The inverse cosine interpolation policy
  typedef LinLogCos InterpPolicy;
};

//! Helper class used to invert the cosine interpolation policy (Lin-LogCos)
template<>
struct InverseInterpPolicy<LinLogCos>
{
  //! The inverse cosine interpolation policy
  typedef LogCosLin InterpPolicy;
};

} // end Utility namespace

//---------------------------------------------------------------------------//
// Template includes
//---------------------------------------------------------------------------//

#include "Utility_CosineInterpolationPolicy_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_COSINE_INTERPOLATION_POLICY_HPP

//---------------------------------------------------------------------------//
// end Utility_CosineInterpolationPolicy.hpp
//---------------------------------------------------------------------------//
