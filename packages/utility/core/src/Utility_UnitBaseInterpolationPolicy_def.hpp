//---------------------------------------------------------------------------//
//!
//! \file   Utility_UnitBaseInterpolationPolicy_def.hpp
//! \author Alex Robinson
//! \brief  Unit base interpolation policy struct definition
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iterator>

// FRENSIE Includes
#include "Utility_InterpolationPolicy.hpp"
#include "Utility_TupleMemberTraits.hpp"
#include "Utility_ContractException.hpp"

namespace Utility{

// Calculate the intermediate range of the Y independent variable
template<typename ZYInterpPolicy, 
	 typename ZXInterpPolicy,
	 typename RXInterpPolicy>
template<typename T>
inline T UnitBaseImpl<ZYInterpPolicy,
		      ZXInterpPolicy,
		      RXInterpPolicy>::calculateIntermediateRange( 
						  const T indep_var_x_0,
						  const T indep_var_x_1,
						  const T indep_var_x,
						  const T indep_var_y_0_range,
						  const T indep_var_y_1_range )
{
  // The interpolation type on the Z variable must be consistent
  testStaticPrecondition( (boost::is_same<typename ZYInterpPolicy::DepVarProcessingTag,typename ZXInterpPolicy::DepVarProcessingTag>::value) );
  // T must be a floating point type
  testStaticPrecondition( (boost::is_floating_point<T>::value) );
  // Make sure the first independent variables are valid
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( indep_var_x_0 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( indep_var_x_1 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( indep_var_x ) );
  testPrecondition( ZXInterpPolicy::isIndepVarInValidRange( indep_var_x_0 ) );
  testPrecondition( indep_var_x_0 < indep_var_x_1 );
  testPrecondition( indep_var_x >= indep_var_x_0 );
  testPrecondition( indep_var_x <= indep_var_x_1 );
  // Make sue the min second independent variables are valid
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( min_indep_var_y_0 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( min_indep_var_y_1 ) );
  testPrecondition( indep_var_y_0_range > 0.0 );
  testPrecondition( indep_var_y_1_range > 0.0 );
  
  return RXInterpPolicy::interpolate( indep_var_x_0,
				      indep_var_x_1,
				      indep_var_x,
				      indep_var_y_0_range,
				      indep_var_y_1_range );
}

// Calculate the intermediate min Y independent variable
template<typename ZYInterpPolicy, 
	 typename ZXInterpPolicy,
	 typename RXInterpPolicy>
template<typename T>
inline T UnitBaseImpl<ZYInterpPolicy,
		      ZXInterpPolicy,
		      RXInterpPolicy>::calculateIntermediateMinYIndepValue( 
						    const T indep_var_x_0,
						    const T indep_var_x_1,
						    const T indep_var_x,
						    const T min_indep_var_y_0,
						    const T min_indep_var_y_1 )
{
  // The interpolation type on the Z variable must be consistent
  testStaticPrecondition( (boost::is_same<typename ZYInterpPolicy::DepVarProcessingTag,typename ZXInterpPolicy::DepVarProcessingTag>::value) );
  // T must be a floating point type
  testStaticPrecondition( (boost::is_floating_point<T>::value) );
  // Make sure the first independent variables are valid
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( indep_var_x_0 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( indep_var_x_1 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( indep_var_x ) );
  testPrecondition( ZXInterpPolicy::isIndepVarInValidRange( indep_var_x_0 ) );
  testPrecondition( indep_var_x_0 < indep_var_x_1 );
  testPrecondition( indep_var_x >= indep_var_x_0 );
  testPrecondition( indep_var_x <= indep_var_x_1 );
  // Make sue the min second independent variables are valid
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( min_indep_var_y_0 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( min_indep_var_y_1 ) );
  testPrecondition( ZYInterpPolicy::isIndepVarInValidRange(min_indep_var_y_0));
  testPrecondition( ZYInterpPolicy::isIndepVarInValidRange(min_indep_var_y_1));
  
  return RXInterpPolicy::interpolate( indep_var_x_0,
				      indep_var_x_1,
				      indep_var_x,
				      min_indep_var_y_0,
				      min_indep_var_y_1 );
}

// Calculate the unit base y independent variable
template<typename ZYInterpPolicy, 
	 typename ZXInterpPolicy,
	 typename RXInterpPolicy>
template<typename T>
inline T UnitBaseImpl<ZYInterpPolicy,
		      ZXInterpPolicy,
		      RXInterpPolicy>::calculateUnitBaseIndepYVariable( 
						    const T indep_var_y,
						    const T min_indep_var_y,
					            const T indep_var_y_range )
{
  // Make sure the independent variables are valid
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( indep_var_y ) );
  testPrecondition( ZYInterpPolicy::isIndepVarInValidRange( indep_var_y ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( min_indep_var_y ) );
  testPrecondition( ZYInterpPolicy::isIndepVarInValidRange( min_indep_var_y ));
  testPrecondition( indep_var_y >= min_indep_var_y );
  // Make sure the range is valid
  testPrecondition( indep_var_y_range > 0.0 );
  testPrecondition( indep_var_y_range >= indep_var_y - min_indep_var_y );
  
  return (indep_var_y - min_indep_var_y)/indep_var_y_range;
}

// Calculate the independent y variable
template<typename ZYInterpPolicy, 
	 typename ZXInterpPolicy,
	 typename RXInterpPolicy>
template<typename T>
inline T UnitBaseImpl<ZYInterpPolicy,
		      ZXInterpPolicy,
		      RXInterpPolicy>::calculateIndepYVariable( 
						 const T unit_base_indep_var_y,
						 const T min_indep_var_y,
						 const T indep_var_y_range )
{
  // Make sure the unit base independent variable is valid
  testPrecondition( unit_base_indep_var_y >= 0.0 );
  testPrecondition( unit_base_indep_var_y <= 1.0 );
  // Make sure the min independent variable is valid
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( min_indep_var_y ) );
  testPrecondition( ZYInterpPolicy::isIndepVarInValidRange( min_indep_var_y ));
  // Make sure the range is valid
  testPrecondition( indep_var_y_range > 0.0 );
  
  return min_indep_var_y + indep_var_y_range*unit_base_indep_var_y;
}

// Interpolate between four points (x0,y00), (x0,y01), (x1,y10), (x0,y11)
template<typename ZYInterpPolicy, 
	 typename ZXInterpPolicy,
	 typename RXInterpPolicy>
template<TupleMember YIndepMember,
	 TupleMember DepMember,
	 typename T, 
	 typename Iterator>
inline T UnitBaseImpl<ZYInterpPolicy,
		      ZXInterpPolicy,
		      RXInterpPolicy>::interpolate(
						  const T indep_var_x_0,
						  const T indep_var_x_1,
						  const T indep_var_x,
						  const T indep_var_y,
						  Iterator start_indep_var_y_0,
						  Iterator end_indep_var_y_0,
						  Iterator start_dep_var_0,
						  Iterator end_dep_var_0,
						  Iterator start_indep_var_y_1,
						  Iterator end_indep_var_y_1,
						  Iterator start_dep_var_1,
						  Iterator start_dep_var_1 )
{
  // The interpolation type on the Z variable must be consistent
  testStaticPrecondition( (boost::is_same<typename ZYInterpPolicy::DepVarProcessingTag,typename ZXInterpPolicy::DepVarProcessingTag>::value) );
  // T must be a floating point type
  testStaticPrecondition( (boost::is_floating_point<T>::value) );
  // The iterator must have T as the value type
  testStaticPrecondition( (boost::is_same<typename TupleMemberTraits<typename std::iterator_traits<Iterator>::value_type,member>::tupleMemberType,T>::value) );
  // Make sure the first independent variables are valid
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( indep_var_x_0 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( indep_var_x_1 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( indep_var_x ) );
  testPrecondition( ZXInterpPolicy::isIndepVarInValidRange( indep_var_x_0 ) );
  testPrecondition( indep_var_x_0 < indep_var_x_1 );
  testPrecondition( indep_var_x >= indep_var_x_0 );
  testPrecondition( indep_var_x <= indep_var_x_1 );
  // Make sure the second independent variables are valid
  testPrecondition( start_indep_var_y_0 != end_indep_var_y_0 );
  testPrecondition( Sort::isSortedAscending<YIndepMember>(start_indep_var_y_0,
							  end_indep_var_y_0) );
  testPrecondition( ZYInterpPolicy::isIndepVarInValidRange( 
				  get<YIndepMember>( start_indep_var_y_0 ) ) );
  testPrecondition( start_indep_var_y_1 != end_indep_var_y_1 );
  testPrecondition( Sort::isSortedAscending<YIndepMember>(start_indep_var_y_1,
							  end_indep_var_y_1) );
  testPrecondition( ZYInterpPolicy::isIndepVarInValidRange( 
				  get<YIndepMember>( start_indep_var_y_1 ) ) );
  // Make sure the dependent variables are valid
  testPrecondition( start_dep_var_y_0 != end_dep_var_y_0 );
  testPrecondition( start_dep_var_y_1 != end_dep_var_y_1 );

  // Calculate the intermediate
  
  // Interpolate on the first grid
  double dep_var_x_0_y = ZYInterpPolicy::interpolate( indep_var_y_0_0,
						      indep_var_y_0_1,
						      indep_var_y,
						      dep_var_0_0,
						      dep_var_0_1 );

  // Interpolate on the second grid
  double dep_var_x_1_y = ZYInterpPolicy::interpolate( indep_var_y_1_0,
						      indep_var_y_1_1,
						      indep_var_y,
						      dep_var_1_0,
						      dep_var_1_1 );

  // Calculate the intermediate range
  double intermediate_range = RXInterpPolicy::interpolate(
							 indep_var_x_0,
							 indep_var_x_1,
							 indep_var_x,
							 indep_var_y_0_range,
						         indep_var_y_1_range );
							  
  // Make sure the intermediate range is valid
  testPostcondition( intermediate_range > 0.0 );
  testPostcondition( intermediate_range >= indep_var_y_0_range );
  testPostcondition( intermediate_range <= indep_var_y_1_range );

  // Conduct the unit base interpolation between grids
  return ZXInterpPolicy::interpolate( indep_var_x_0,
				      indep_var_x_1,
				      indep_var_x,
				      dep_var_x_0_y*indep_var_y_0_range,
				      dep_var_x_1_y*indep_var_y_1_range )/
    intermediate_range;
}

// Interpolate between four processed points
template<typename ZYInterpPolicy, 
	 typename ZXInterpPolicy,
	 typename RXInterpPolicy>
template<typename T>
inline T UnitBaseImpl<ZYInterpPolicy,ZXInterpPolicy>::::interpolate( 
					 const T processed_indep_var_x_0,
					 const T processed_indep_var_x_1,
					 const T processed_indep_var_x,
					 const T processed_indep_var_y_0_0,
					 const T indep_var_y_0_range,
					 const T processed_dep_var_0_0,
					 const T processed_slope_0,
					 const T processed_indep_var_y_1_0,
					 const T indep_var_y_1_range,
					 const T processed_dep_var_1_0,
					 const T processed_slope_1,
					 const T processed_indep_var_y )
{
  // The interpolation type on the Z variable must be consistent
  testStaticPrecondition( (boost::is_same<typename ZYInterpPolicy::DepVarProcessingTag,typename ZXInterpPolicy::DepVarProcessingTag>::value) );
  // T must be a floating point type
  testStaticPrecondition( (boost::is_floating_point<T>::value) );
  // Make sure the first processed independent variables are valid
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf(
						    processed_indep_var_x_0) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf(
						    processed_indep_var_x_1) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf(
						      processed_indep_var_x) );
  testPrecondition( processed_indep_var_x_0 < processed_indep_var_x_1 );
  testPrecondition( processed_indep_var_x >= processed_indep_var_x_0 );
  testPrecondition( processed_indep_var_x <= processed_indep_var_x_1 );
  // Make sure the second independent variables are valid
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf(
						 processed_indep_var_y_0_0 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( 
						 processed_indep_var_y_1_0 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( 
						     processed_indep_var_y ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( processed_slope_0 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf( processed_slope_1 ) );
  testPrecondition( processed_indep_var_y >= indep_var_y_0_0 );
  testPrecondition( processed_indep_var_y >= indep_var_y_1_0 );
  // Make sure the second independent variable ranges are valid
  testPrecondition( indep_var_y_0_range > 0.0 );
  testPrecondition( indep_var_y_1_range > 0.0 );
  // Make sure the dependent variables are valid
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf(
						     processed_dep_var_0_0 ) );
  testPrecondition( !Teuchos::ScalarTraits<T>::isnaninf(
						     processed_dep_var_1_0 ) );
  
  // Interpolate on the first grid
  double dep_var_x_0_y = ZYInterpPolicy::interpolate(processed_indep_var_y_0_0,
						     processed_indep_var_y,
						     processed_dep_var_0_0,
						     processed_slope_0 );

  // Process the value on the first grid
  double processed_dep_var_x_0_y = 
    ZYInterpPolicy::processDepVar( dep_var_x_0_y*indep_var_y_0_range );

  // Interpolate on the second grid
  double dep_var_x_1_y = ZYInterpPolicy::interpolate(processed_indep_var_y_1_0,
						     processed_indep_var_y,
						     processed_dep_var_1_0,
						     processed_slope_1 );

  // Process the value on the second grid
  double processed_dep_var_x_1_y = 
    ZYInterpPolicy::processDepVar( dep_var_x_1_y*indep_var_y_1_range );
  
  // Determine the processed slope between the grids
  double processed_slope = (processed_dep_var_x_1_y - processed_dep_var_x_0_y)/
    (processed_indep_var_x_1 - processed_indep_var_x_0);

  // Calculate the intermediate range
  double processed_indep_var_y_0_range = 
    ZXInterpPolicy::processDepVar( indep_var_y_0_range );

  double processed_indep_var_y_1_range = 
    ZXInterpPolicy::processDepVar( indep_var_y_1_range );
  
  double processed_range_slope = (processed_indep_var_y_1_range - 
				  processed_indep_var_y_0_range )/
    (processed_indep_var_x_1 - processed_indep_var_x_0);
  
  double intermediate_range = RXInterpPolicy::interpolate(
						 processed_indep_var_x_0,
						 processed_indep_var_x,
						 processed_indep_var_y_0_range,
						 processed_range_slope );

  // Make sure the intermediate range is valid
  testPostcondition( intermediate_range > 0.0 );
  testPostcondition( intermediate_range >= indep_var_y_0_range );
  testPostcondition( intermediate_range <= indep_var_y_1_range );

  // Conduct the unit base interpolation between grids
  return ZXInterpPolicy::interpolate( processed_indep_var_x_0,
				      processed_indep_var_x,
				      processed_dep_var_x_0_y,
				      processed_slope )/intermediate_range;
}

} // end Utility namespace

//---------------------------------------------------------------------------//
// end Utility_UnitBaseInterpolationPolicy_def.hpp
//---------------------------------------------------------------------------//
