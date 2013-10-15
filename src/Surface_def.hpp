//---------------------------------------------------------------------------//
//!
//! \file   Surface_def.hpp
//! \author Alex Robinson
//! \brief  Surface class definition
//!
//---------------------------------------------------------------------------//

#ifndef SURFACE_DEF_HPP
#define SURFACE_DEF_HPP

// Std Lib Includes
#include <sstream>

// Trilinos Includes
#include <Teuchos_SerialDenseHelpers.hpp>

// FACEMC Includes
#include "ContractException.hpp"
#include "VectorHelpers.hpp"
#include "MatrixHelpers.hpp"
#include "LinearAlgebraAlgorithms.hpp"

namespace FACEMC{

// General Surface Constructor
template<typename OrdinalType, typename ScalarType>
Surface<OrdinalType,ScalarType>::Surface( OrdinalType id,
					  ScalarType a,
					  ScalarType b,
					  ScalarType c,
					  ScalarType d,
					  ScalarType e,
					  ScalarType f,
					  ScalarType g,
					  ScalarType h,
					  ScalarType j,
					  ScalarType k,
					  ScalarType tolerance_ratio )
  : PrintableObject( "ThreeSpaceSurface" ),
    ThreeSpaceObject( THREE_SPACE_SURFACE ),
    d_id( id ),
    d_definition( Teuchos::tuple( a, b, c, d, e, f, g, h, j, k ) ),
    d_tolerance( ST::zero() ),
    d_symmetric( false ),
    d_planar( false )
{ 
  // Make sure the id is valid
  testPrecondition( id <= OT::max() && id >= OT::zero() )
  // Make sure the correct constructor has been used
  testPrecondition( (a != ST::zero() || b != ST::zero() || c != ST::zero()) && 
		    (d != ST::zero() || e != ST::zero() || f != ST::zero()) );

  // Determine floating point tolerance
  setTolerance( tolerance_ratio );
}

// Symmetric 2nd order surface constructor
template<typename OrdinalType, typename ScalarType>
Surface<OrdinalType,ScalarType>::Surface( OrdinalType id,
					  ScalarType a,
					  ScalarType b,
					  ScalarType c,
					  ScalarType g,
					  ScalarType h,
					  ScalarType j,
					  ScalarType k,
					  ScalarType tolerance_ratio )
  : PrintableObject( "ThreeSpaceSurface" ),
    ThreeSpaceObject( THREE_SPACE_SURFACE ),
    d_id( id ),
    d_definition( Teuchos::tuple( a, b, c, ST::zero(), ST::zero(), ST::zero(), g, h, j, k ) ),
    d_tolerance( ST::zero() ),
    d_symmetric( true ),
    d_planar( false )
{
  // Make sure the id is valid
  testPrecondition( id <= OT::max() && id >= OT::zero() );
  // Make sure the correct constructor has been used
  testPrecondition( a != ST::zero() || b != ST::zero() || c != ST::zero() );

  // Determine floating point tolerance
  setTolerance( tolerance_ratio );
}

// Planar Surface constructor
template<typename OrdinalType, typename ScalarType>
Surface<OrdinalType,ScalarType>::Surface( OrdinalType id,
					  ScalarType g,
					  ScalarType h,
					  ScalarType j,
					  ScalarType k,
					  ScalarType tolerance_ratio )
  : PrintableObject( "ThreeSpaceSurface" ),
    ThreeSpaceObject( THREE_SPACE_SURFACE ),
    d_id( id ),
    d_definition( Teuchos::tuple( ST::zero(), ST::zero(), ST::zero(), ST::zero(), ST::zero(), ST::zero(), g, h, j, k ) ),
    d_tolerance( ST::zero() ),
    d_symmetric( true ),
    d_planar( true )
{
  // Make sure that the id is valid
  testPrecondition( id <= OT::max() && id >= OT::zero() );
  // Make sure that the surface is valid
  testPrecondition( g != ST::zero() || h != ST::zero() || j != ST::zero() );

  // Scale the surface so that the linear term vector is a unit vector (this
  // will not change the nature of the surface - it will merely make
  // transformations using this surface cleaner).
  ScalarType scale_factor = getLinearTermVector().normTwo();
  d_definition[6] /= scale_factor;
  d_definition[7] /= scale_factor;
  d_definition[8] /= scale_factor;
  d_definition[9] /= scale_factor;

  // Determine floating point tolerance
  setTolerance( tolerance_ratio );
}

// Construct surface by translating another surface
template<typename OrdinalType, typename ScalarType>
Surface<OrdinalType,ScalarType>::Surface( 
        OrdinalType id,
	const Surface<OrdinalType,ScalarType> &original_surface,
	const Vector<ScalarType> &translation_vector )
  : PrintableObject( "ThreeSpaceSurface" ),
    ThreeSpaceObject( THREE_SPACE_SURFACE ),
    d_id( id ),
    d_definition( original_surface.d_definition ),
    d_tolerance( original_surface.d_tolerance ),
    d_symmetric( original_surface.d_symmetric ),
    d_planar( original_surface.d_planar )
{
  // Make sure that the id is valid
  testPrecondition( id <= OT::max() && id >= OT::zero() );
  // Make sure that the translation vector is valid
  testPrecondition( translation_vector.normTwo() > ST::zero() );

  // Matrix form of surface eqn: x^tAx+b^tx+k => ^t indicates the transpose
  // Translate the surface using the given translation vector:
  //  A' = A
  //  b' = (A + A^t)x0 => x0 is translation vector
  //  k' = x0^tAx0 + b^tx0 + k
  Vector<ScalarType> b = getLinearTermVector();
  
  if( d_planar )
  {
    // k' = b^tx0 + k
    d_definition[9] += translation_vector.dot( b );
  }
  else // 2nd order surface
  {
    Matrix<ScalarType> A = getQuadraticFormMatrix();
    Vector<ScalarType> Ax0, Atx0;

    Ax0.multiply( 1.0, A, false, translation_vector, 0.0 );
    
    Atx0.multiply( 1.0, A, true, translation_vector, 0.0 );
    
    // g'
    d_definition[6] = Ax0[0] + Atx0[0];
    // h'
    d_definition[7] = Ax0[1] + Atx0[1];
    // j'
    d_definition[8] = Ax0[2] + Atx0[2];
    // k'
    d_definition[9] += translation_vector.dot( Ax0 ) + 
      translation_vector.dot( b );
  }

  // The new surface must be distinct from the original
  testPostcondition( (getQuadraticFormMatrix() != 
		      original_surface.getQuadraticFormMatrix()) ||
		     (getLinearTermVector() !=
		      original_surface.getLinearTermVector()) ||
		     (getConstantTerm() !=
		      original_surface.getConstantTerm()) );
}

// Construct surface by rotating another surface
template<typename OrdinalType, typename ScalarType>
Surface<OrdinalType,ScalarType>::Surface( 
	   OrdinalType id,
	   const Surface<OrdinalType,ScalarType> &original_surface,
	   const Matrix<ScalarType> &rotation_matrix )
  : PrintableObject( "ThreeSpaceSurface" ),
    ThreeSpaceObject( THREE_SPACE_SURFACE ),
    d_id( id ),
    d_definition( original_surface.d_definition ),
    d_tolerance( original_surface.d_tolerance ),
    d_symmetric( original_surface.d_symmetric ),
    d_planar( original_surface.d_planar )
{
  // Make sure that the id is valid
  testPrecondition( id <= OT::max() && id >= OT::zero() );
  // Make sure that the rotation matrix is valid (3x3 and orthonormal)
  testPrecondition( rotation_matrix.isOrthonormal() );

  remember( int multiply_success = 0 ); // successful multiplication returns 0

  // Matrix form of surface eqn: x^tAx+b^tx+k => ^t indicates the transpose
  // Translate the surface using the given rotation matrix:
  //  A' = R^tAR => R is the rotation matrix
  //  b' = R^tb
  //  k' = k
  Vector<ScalarType> b = getLinearTermVector(), b_prime;
  
  b_prime.multiply( 1.0, rotation_matrix, true, b, 0.0 );

  // g'
  d_definition[6] = b_prime[0];
  // h'
  d_definition[7] = b_prime[1];
  // j'
  d_definition[8] = b_prime[2];

  if( !d_planar ) // second order surface
  {
    Matrix<ScalarType> A = getQuadraticFormMatrix(), A_prime, AR;
    
    AR.multiply( 1.0, A, false, rotation_matrix, false, 0.0 );
    A_prime.multiply( 1.0, rotation_matrix, true, AR, true, 0.0 );
  
    // a'
    d_definition[0] = A_prime( 0, 0 );
    // b'
    d_definition[1] = A_prime( 1, 1 );
    // c'
    d_definition[2] = A_prime( 2, 2 );
    // d'
    d_definition[3] = A_prime( 0, 1 )*2;
    // e'
    d_definition[4] = A_prime( 1, 2 )*2;
    // f'
    d_definition[5] = A_prime( 0, 2 )*2;

    // Eliminate small values (~0.0)
    filterDefinition();

    // check if the surface symmetry has changed
    checkSymmetry();
  }

  // The new surface must be distinct from the original
  testPostcondition( (getQuadraticFormMatrix() != 
		      original_surface.getQuadraticFormMatrix()) ||
		     (getLinearTermVector() !=
		      original_surface.getLinearTermVector()) ||
		     (getConstantTerm() !=
		      original_surface.getConstantTerm()) );
}

// Construct surface by conducting a general transform on anothe surface
template<typename OrdinalType, typename ScalarType>
Surface<OrdinalType,ScalarType>::Surface( 
	OrdinalType id,
	const Surface<OrdinalType,ScalarType> &original_surface,
	const Matrix<ScalarType> &rotation_matrix,
        const Vector<ScalarType> &translation_vector )
  : PrintableObject( "ThreeSpaceSurface" ),
    ThreeSpaceObject( THREE_SPACE_SURFACE ),
    d_id( id ),
    d_definition( original_surface.d_definition ),
    d_tolerance( original_surface.d_tolerance ),
    d_symmetric( original_surface.d_symmetric ),
    d_planar( original_surface.d_planar )
{
  // Make sure that the id is valid
  testPrecondition( id <= OT::max() && id >= OT::zero() );
  // Make sure that the rotation matrix is valid (3x3 and orthonormal)
  testPrecondition( rotation_matrix.isOrthonormal() );
  // Make sure that the translation vector is valid
  testPrecondition( translation_vector.normTwo() > ST::zero() );

  // Matrix form of surface eqn: x^tAx+b^tx+k => ^t indicates the transpose
  // Translate the surface using the given rot. matrix & trans. vector:
  //  A' = R^tAR => R is the rotation matrix
  //  b' = R^t(A + A^t)x0 + R^tb => x0 is the translation vector
  //  k' = x0^tAx0 + b^tx0 + k
  Vector<ScalarType> b = getLinearTermVector(), b_prime;

  b_prime.multiply( 1.0, rotation_matrix, true, b, 0.0 );

  if( d_planar )
  {
    // g'
    d_definition[6] = b_prime[0];
    // h'
    d_definition[7] = b_prime[1];
    // j'
    d_definition[8] = b_prime[2];
    // k'
    d_definition[9] += translation_vector.dot( b );
  }
  else // second order surface
  {
    Matrix<ScalarType> A = getQuadraticFormMatrix(), A_prime;
    Vector<ScalarType> AR, Ax0, Atx0;
    
    Ax0.multiply( 1.0, A, false, translation_vector, 0.0 );
    
    Atx0.multiply( 1.0, A, true, translation_vector, 0.0 );
    
    // Add extra terms to b'
    b_prime.multiply( 1.0, rotation_matrix, true, Ax0, 1.0 );
    b_prime.multiply( 1.0, rotation_matrix, true, Atx0, 1.0 );
    
    // Compute new quadratic form matrix
    AR.multiply( 1.0, A, false, rotation_matrix, false, 0.0 );
    A_prime.multiply( 1.0, rotation_matrix, true, AR, false, 0.0 );
  
    // a'
    d_definition[0] = A_prime( 0, 0 );
    // b'
    d_definition[1] = A_prime( 1, 1 );
    // c'
    d_definition[2] = A_prime( 2, 2 );
    // d'
    d_definition[3] = A_prime( 0, 1 )*2;
    // e'
    d_definition[4] = A_prime( 1, 2 )*2;
    // f'
    d_definition[5] = A_prime( 0, 2 )*2;
    // g'
    d_definition[6] = b_prime[0];
    // h'
    d_definition[7] = b_prime[1];
    // j'
    d_definition[8] = b_prime[2];
    // k'
    d_definition[9] += translation_vector.dot( Ax0 ) + 
      translation_vector.dot( b );

    // Eliminate small values (~0.0)
    filterDefinition();

    // check if the surface symmetry has changed
    checkSymmetry();
  }
  
  // The new surface must be distinct from the original
  testPostcondition( (getQuadraticFormMatrix() != 
		      original_surface.getQuadraticFormMatrix()) ||
		     (getLinearTermVector() !=
		      original_surface.getLinearTermVector()) ||
		     (getConstantTerm() !=
		      original_surface.getConstantTerm()) );
}

// Return if the point is on the surface
template<typename OrdinalType, typename ScalarType>
bool Surface<OrdinalType,ScalarType>::isOn( const ScalarType x, 
					    const ScalarType y, 
					    const ScalarType z ) const
{
  ScalarType surface_evaluation = evaluateSurface( x, y, z );

  return ST::magnitude( surface_evaluation ) < d_tolerance;
}

// Return if the point is on the surface
template<typename OrdinalType, typename ScalarType>
bool Surface<OrdinalType,ScalarType>::isOn( 
					const Vector<ScalarType> &point ) const
{
  return isOn( point[0], point[1], point[2] );
}

// Return if the surface is planar
template<typename OrdinalType, typename ScalarType>
bool Surface<OrdinalType,ScalarType>::isPlanar() const
{
  return d_planar;
}

// Return if the surface is symmetric
template<typename OrdinalType, typename ScalarType>
bool Surface<OrdinalType,ScalarType>::isSymmetric() const
{
  return d_symmetric;
}

// Return the surface id
template<typename OrdinalType, typename ScalarType>
OrdinalType Surface<OrdinalType,ScalarType>::getId() const
{
  return d_id;
}

// Return the sense of a point with respect to the surface
template<typename OrdinalType, typename ScalarType>
SurfaceSense Surface<OrdinalType,ScalarType>::getSenseOfPoint( 
						     const ScalarType x, 
						     const ScalarType y, 
						     const ScalarType z ) const
{
  ScalarType surface_evaluation = evaluateSurface( x, y, z );
  
  if( surface_evaluation > d_tolerance )
    return POS_SURFACE_SENSE;
  else if( surface_evaluation < -d_tolerance )
    return NEG_SURFACE_SENSE;
  else 
    return ON_SURFACE;
}

// Return the sense of a point with respect to the surface
template<typename OrdinalType, typename ScalarType>
SurfaceSense Surface<OrdinalType,ScalarType>::getSenseOfPoint( 
				        const Vector<ScalarType> &point ) const
{
  return getSenseOfPoint( point[0], point[1], point[2] );
}
    
// Return the unit normal from the surface at a point on the surface,
// pointing in the direction of the desired sense
template<typename OrdinalType, typename ScalarType>
Vector<ScalarType> Surface<OrdinalType,ScalarType>::getUnitNormalAtPoint( 
					       const ScalarType x, 
					       const ScalarType y, 
					       const ScalarType z,
					       const SurfaceSense sense ) const
{
  // The point must be on the surface
  testPrecondition( isOn( x, y, z ) );
  // The surface sense must be either pos or neg (not on)
  testPrecondition( sense != ON_SURFACE );

  Vector<ScalarType> unit_normal;
  
  if( d_planar )
  {
    // dS/dx: g
    unit_normal[0] = d_definition[6];
    // dS/dy: h
    unit_normal[1] = d_definition[7];
    // dS/dz: j
    unit_normal[2] = d_definition[8];
    // Note: planes are normalized at construction
  }
  else
  {
    // dS/dx: 2ax+dy+fz+g
    unit_normal[0] = 2*d_definition[0]*x + d_definition[3]*y + 
      d_definition[5]*z + d_definition[6];
    
    // dS/dy: 2by+dx+ez+h
    unit_normal[1] = 2*d_definition[1]*y + d_definition[3]*x + 
      d_definition[4]*z + d_definition[7];
    
    // dS/dz: 2cz+ey+fx+j
    unit_normal[2] = 2*d_definition[2]*z + d_definition[4]*y + 
      d_definition[5]*x + d_definition[8];
    
    // Normalize the vector
    unit_normal.normalize();
  }
  
  // Point the vector in the direction of the desired half-space
  SurfaceSense sense_of_norm_tip = getSenseOfPoint( x + unit_normal[0],
						    y + unit_normal[1],
						    z + unit_normal[2] );

  // Potentially reverse the direction of the unit normal
  if( sense_of_norm_tip != sense )
    unit_normal *=-1.0;

  // The sense of the normal vector tip must not be ON_SURFACE
  testPostcondition( sense_of_norm_tip != ON_SURFACE );
  // The unit normal calculated must be valid
  testPostcondition( ST::magnitude( unit_normal.normTwo() - 1.0 ) < 
		     ST::prec() );

  return unit_normal;
}

// Return the unit normal from the surface at a point on the surface,
// pointing in the direction of the desired sense
template<typename OrdinalType, typename ScalarType>
Vector<ScalarType> Surface<OrdinalType,ScalarType>::getUnitNormalAtPoint( 
					       const Vector<ScalarType> &point,
					       const SurfaceSense sense ) const
{
  return getUnitNormalAtPoint( point[0], point[1], point[2] );
}

// Return the quadratic form matrix of the surface
template<typename OrdinalType, typename ScalarType>
Matrix<ScalarType>
Surface<OrdinalType,ScalarType>::getQuadraticFormMatrix() const
{
  return Matrix<ScalarType>( 
		       d_definition[0], 
		       d_definition[3]/2, d_definition[1],
		       d_definition[5]/2, d_definition[4]/2, d_definition[2] );
}

// Return the linear term vector of the surface
template<typename OrdinalType, typename ScalarType>
Vector<ScalarType> Surface<OrdinalType,ScalarType>::getLinearTermVector() const
{
  return Vector<ScalarType>( d_definition[6], 
			     d_definition[7], 
			     d_definition[8] );
}

//! Return the constant term of the surface
template<typename OrdinalType, typename ScalarType>
ScalarType Surface<OrdinalType,ScalarType>::getConstantTerm() const
{
  return d_definition[9];
}

// Print method that defines the behavior of the std::stream << operator
template<typename OrdinalType, typename ScalarType>
void Surface<OrdinalType,ScalarType>::print( std::ostream &os ) const
{
  std::stringstream ss;
  if( ST::magnitude( d_definition[0] ) > ST::prec() )
    ss << d_definition[0] << "x^2 ";
  if( ST::magnitude( d_definition[1] ) > ST::prec() )
    ss << d_definition[1] << "y^2 ";
  if( ST::magnitude( d_definition[2] ) > ST::prec() )
    ss << d_definition[2] << "z^2 ";
  if( ST::magnitude( d_definition[3] ) > ST::prec() )
    ss << d_definition[3] << "xy ";
  if( ST::magnitude( d_definition[4] ) > ST::prec() )
    ss << d_definition[4] << "yz ";
  if( ST::magnitude( d_definition[5] ) > ST::prec() )
    ss << d_definition[5] << "xz ";
  if( ST::magnitude( d_definition[6] ) > ST::prec() )
    ss << d_definition[6] << "x ";
  if( ST::magnitude( d_definition[7] ) > ST::prec() )
    ss << d_definition[7] << "y ";
  if( ST::magnitude( d_definition[8] ) > ST::prec() )
    ss << d_definition[8] << "z ";
  if( ST::magnitude( d_definition[9] ) > ST::prec() )
    ss << d_definition[9];
  
  std::string output = ss.str();
  
  // remove trailing spaces
  if( output[output.size()-1] == ' ' )
    output.erase( output.size()-1, 1 );
  
  // replace all remaining spaces with '+'
  unsigned white_space_loc = output.find( " " );
  
  while( white_space_loc < output.size() )
  {
    output[white_space_loc] = '+';
    
    white_space_loc = output.find( " ", white_space_loc );
  }

  // print the string
  os << "Id: " << d_id << std::endl;
  os << "Definition: " << output << " = 0 " << std::endl;  
}

// Calculate tolerance based on tolerance ratio
template<typename OrdinalType, typename ScalarType>
void Surface<OrdinalType,ScalarType>::setTolerance( ScalarType tolerance_ratio)
{
  typename ST::magnitudeType max = ST::zero();
  
  for( tupleIndex i = Tuple_OT::zero(); i < 10; ++i )
  {
    if( max < ST::magnitude( d_definition[i] ) )
      max = ST::magnitude( d_definition[i] );
  
    d_tolerance = tolerance_ratio*max;
  }
}

// Filter out small values that can be calculated from a transform
template<typename OrdinalType, typename ScalarType>
void Surface<OrdinalType,ScalarType>::filterDefinition(ScalarType cutoff_value)
{
  for( tupleIndex i = Tuple_OT::zero(); i < 10; ++i )
  {
    if( ST::magnitude( d_definition[i] ) < cutoff_value )
      d_definition[i] = ST::zero();
  }
}

// Check if the surface symmetry has changed after a transform
template<typename OrdinalType, typename ScalarType>
void Surface<OrdinalType,ScalarType>::checkSymmetry()
{
  if( d_symmetric )
  {
    if( (d_definition[3] != ST::zero()) || (d_definition[4] != ST::zero()) ||
	(d_definition[5] != ST::zero()) )
      d_symmetric = false; 
  }
  else
  {
    if( (d_definition[3] == ST::zero()) && (d_definition[4] == ST::zero()) &&
	(d_definition[5] == ST::zero()) )
      d_symmetric = true; 
  }
}

// Evaluate the surface definition at a point
template<typename OrdinalType, typename ScalarType>
ScalarType Surface<OrdinalType,ScalarType>::evaluateSurface( 
						     const ScalarType x, 
						     const ScalarType y, 
						     const ScalarType z ) const
{
  if( d_planar )
    return evaluatePlanarSurface( x, y, z );
  else if( d_symmetric )
    return evaluateSymmetricSecondOrderSurface( x, y, z );
  else
    return evaluateGeneralSecondOrderSurface( x, y, z );
}

// Evaluate planar surface
template<typename OrdinalType, typename ScalarType>
ScalarType Surface<OrdinalType,ScalarType>::evaluatePlanarSurface( 
						     const ScalarType x, 
						     const ScalarType y, 
						     const ScalarType z ) const
{
  return d_definition[6]*x + d_definition[7]*y + d_definition[8]*z + 
    d_definition[9];
}

template<typename OrdinalType, typename ScalarType>
ScalarType 
Surface<OrdinalType,ScalarType>::evaluateSymmetricSecondOrderSurface( 
						     const ScalarType x,
						     const ScalarType y,
						     const ScalarType z ) const
{
  return d_definition[0]*x*x + d_definition[1]*y*y + d_definition[2]*z*z +
    d_definition[6]*x + d_definition[7]*y + d_definition[8]*z + 
    d_definition[9];
}

template<typename OrdinalType, typename ScalarType>
ScalarType 
Surface<OrdinalType,ScalarType>::evaluateGeneralSecondOrderSurface( 
						     const ScalarType x, 
						     const ScalarType y, 
						     const ScalarType z ) const
{
  return d_definition[0]*x*x + d_definition[1]*y*y + d_definition[2]*z*z +
    d_definition[3]*x*y + d_definition[4]*y*z + d_definition[5]*x*z +
    d_definition[6]*x + d_definition[7]*y + d_definition[8]*z + 
    d_definition[9];
}

} // end FACEMC namespace

#endif // end SURFACE_DEF_HPP

//---------------------------------------------------------------------------//
// end Surface_def.hpp
//---------------------------------------------------------------------------//

