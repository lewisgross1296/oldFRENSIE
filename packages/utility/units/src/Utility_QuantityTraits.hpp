//---------------------------------------------------------------------------//
//!
//! \file   Utility_QuantityTraits.hpp
//! \author Alex Robinson
//! \brief  Quantity traits class specializations
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_QUANTITY_TRAITS_HPP
#define UTILITY_QUANTITY_TRAITS_HPP

// Boost Includes
#include <boost/units/quantity.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/static_rational.hpp>
#include <boost/units/limits.hpp>
#include <boost/units/io.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

// FRENSIE Includes
#include "Utility_QuantityTraitsDecl.hpp"

namespace Utility{

/*! \brief The QuantityTraitsFloatingHelper partial specialization for 
 * boost::units::quantity with floating point types.
 * \ingroup quantity_traits
 */
template<typename Unit, typename T>
struct QuantityTraitsHelper<boost::units::quantity<Unit,T>,typename boost::enable_if<boost::is_floating_point<T> >::type>
{
private:
  typedef boost::units::quantity<Unit,T> QuantityType;
  typedef T RawType;

public:
  static inline QuantityType epsilon()
  { return QuantityType::from_value( std::numeric_limits<RayType>::epsilon() ); }
  
  static inline QuantityType inf()
  { return QuantityType::from_value( std::numeric_limits<RawType>::infinity() ); }

  static inline QuantityType nan()
  { return QuantityType::from_value( std::numeric_limits<RawType>::quiet_NaN() ); }

  static inline bool isnaninf( const QuantityType& a )
  { return Teuchos::ScalarTraits<RawType>::isnaninf( a.value() ); }
};

/*! \brief The partial specialization of QuantityTraits for arithmetic 
 * boost::units::quantity
 * \ingroup quantity_traits
 */
template<typename Unit, typename T>
struct QuantityTraits<boost::units::quantity<Unit,T>, typename boost::enable_if<boost::is_arithmetic<T> >::type> : public QuantityTraitsFloatingHelper<boost::units::quantity<Unit,T> >
{
  typedef Unit UnitType;
  typedef T RawType;
  typedef boost::units::quantity<Unit,T> QuantityType;
  typedef boost::is_floating_point<T> is_floating_point;

  template<boost::units::integer_type N, boost::units::integer_type D = 1>
  struct GetQuantityToPowerType
  {
    typedef typename boost::units::power_typeof_helper<QuantityType,boost::units::static_rational<N,D> >::type type;
  };

  static inline QuantityType zero()
  { return QuantityType::from_value( RawType(0) ); }

  static inline QuantityType one()
  { return QuantityType::from_value( RawType(1) ); }

  //! Possible bug in boost::units::sqrt
  static inline typename GetQuantityToPowerType<1,2>::type sqrt( const QuantityType& quantity )
  { 
    return GetQuantityToPowerType<1,2>::type::from_value( std::sqrt( quantity.value() ) ); 
    // return boost::units::sqrt( quantity )
  }

  template<boost::units::integer_type N, boost::units::integer_type D>
  static inline typename GetQuantityToPowerType<N,D>::type rpow( const QuantityType& quantity )
  { return boost::units::pow<boost::units::static_rational<N,D> >( quantity ); }
  //! Potentially dangerous to initialize quantities in this way!
  static inline QuantityType initializeQuantity( const RawType& raw_quantity )
  { return QuantityType::from_value( raw_quantity ); }

  //! Get the raw value of a quantity
  static inline const RawType& getRawQuantity( const QuantityType& quantity )
  { return quantity.value(); }

  //! Set the value of a quantity (potentially dangerous!)
  static inline void setQuantity( QuantityType& quantity,
				  const RawType& raw_quantity )
  { quantity = QuantityType::from_value( raw_quantity ); }
};

/*! \brief The QuantityTraitsFloatingHelper partial specialization for floating
 * point types (no units).
 * \ingroup quantity_traits
 */
template<typename T>
struct QuantityTraitsFloatingHelper<T,typename boost::enable_if<boost::is_floating_point<T> >::type>
{
private:
  typedef T QuantityType;

public:

  static inline QuantityType epsilon()
  { return std::numeric_limits<QuantityType>::epsilon() ); }

  static inline QuantityType inf()
  { return std::numeric_limits<QuantityType>::infinity(); }

  static inline QuantityType nan()
  { return std::numeric_limits<QuantityType>::quiet_NaN(); }

  
  static inline bool isnaninf( const QuantityType& a )
  {
    // Check if NaN
    if( isnanIEEE( a ) )
      return true;

    // Check if Inf
    else if( isInfIEEE( a ) )
      return true;

    // Check if Inf - basic
    else if( isInf( a ) )
      return true;

    // Not NaN or Inf
    else
      return false;
  }

protected:

  /*! Test for quiet NaN using IEEE standard 
   * \details Only guaranteed to work with compliant compilers.
   */
  static inline bool isNaNIEEE( const QuantityType& a )
  { 
    const QuantityType tol = 0.0;
    
    if( !(x <= tol) && !(x > tol) )
      return true;
    else
      return false;
  }

  /*! Test for Inf using IEEE standard
   * \details Only guaranteed to work with compliant compilers.
   */
  static inline bool isInfIEEE( const QuantityType& a )
  {
    // Inf*0.0 -> NaN
    const QuantityType possible_nan = static_cast<QuantityType>(0.0)*a;

    return isNaNIEEE( possible_nan );
  }

  /*! Test for Inf using comparison
   * \details This is not the IEEE recommended test
   */
  static inline isInf( const QuantityType& a )
  {
    if( a == std::numeric_limits<QuantityType>::infinity() ||
        a == -std::numeric_limits<QuantityType>::infinity() )
      return true;
    else
      return false;
  }
};

/*! The specialization of QuantityTraits for all arithmetic types (no units).
 *
 * Note that having no units is different from a dimensionless unit.
 * \ingroup quantity_traits
 */
template<typename T>
struct QuantityTraits<T,typename boost::enable_if<boost::is_arithmetic<T> >::type> : public QuantityTraitsFloatingHelper<T>
{
  typedef void Unit;
  typedef T RawType;
  typedef T QuantityType;
  typedef boost::is_floating_point<T> is_floating_point;

  template<boost::units::integer_type N, boost::units::integer_type D = 1>
  struct GetQuantityToPowerType
  { typedef QuantityType type; };

  static inline QuantityType zero()
  { return RawType(0); }

  static inline QuantityType one()
  { return RawType(1); }

  static inline QuantityType sqrt( const QuantityType quantity )
  { return std::sqrt( quantity ); }

  template<boost::units::integer_type N, boost::units::integer_type D>
  static inline QuantityType rpow( const QuantityType quantity )
  { return boost::units::pow<boost::units::static_rational<N,D> >( quantity ); }
  
  static inline QuantityType initializeQuantity( const RawType& raw_quantity )
  { return raw_quantity; }

  static inline const RawType& getRawQuantity( const QuantityType& quantity )
  { return quantity; }

  static inline void setQuantity( QuantityType& quantity,
				  const RawType& raw_quantity )
  { quantity = raw_quantity; }
};

} // end Utility namespace

#endif // end UTILITY_QUANTITY_TRAITS_HPP

//---------------------------------------------------------------------------//
// end Utility_QuantityTraits.hpp
//---------------------------------------------------------------------------//
