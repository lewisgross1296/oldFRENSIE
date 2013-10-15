//---------------------------------------------------------------------------//
//!
//! \file   tstSurface.cpp
//! \author Alex Robinson
//! \brief  Surface class unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_RCP.hpp>

// FACEMC Includes
#include "FACEMC_UnitTestHarnessExtensions.hpp"
#include "Surface.hpp"
#include "Vector.hpp"
#include "VectorHelpers.hpp"
#include "Matrix.hpp"
#include "MatrixHelpers.hpp"

//---------------------------------------------------------------------------//
// Instantiation Macros.
//---------------------------------------------------------------------------//
#define UNIT_TEST_INSTANTIATION( type, name ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( type, name, short, float ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( type, name, short, double ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( type, name, int, float ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( type, name, int, double ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( type, name, long, float ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( type, name, long, double ) 

//---------------------------------------------------------------------------//
// Testing Info.
//---------------------------------------------------------------------------//
#define TOL 1e-12

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that a surface can determine if a point is on it
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Surface, 
				   isOn,
				   OrdinalType,
				   ScalarType )
{
  // Spherical surface centered at (-10,-10,-10) with radius 10
  FACEMC::Surface<OrdinalType,ScalarType> sphere( 0,
						  1, 1, 1,
						  20, 20, 20,
						  10*10+10*10+10*10 - 10*10 );

  FACEMC::Vector<ScalarType> point_1( 0.0, -10.0, -10.0 );
  FACEMC::Vector<ScalarType> point_2( -10.0, -10.0, -10.0 ); 

  TEST_ASSERT( sphere.isOn( point_1 ) );
  TEST_ASSERT( !sphere.isOn( point_2 ) );
}

UNIT_TEST_INSTANTIATION( Surface, isOn );

//---------------------------------------------------------------------------//
// Check that a surface returns if it is planar or not
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Surface, 
				   isPlanar,
				   OrdinalType,
				   ScalarType )
{
  // Generic surface
  FACEMC::Surface<OrdinalType,ScalarType> surface( 0, 
						   1, 2, 3,
						   4, 5, 6,
						   7, 8, 9,
						   10 );
  
  TEST_ASSERT( !surface.isPlanar() );

  // Spherical surface centered at (-10,-10,-10) with radius 10
  FACEMC::Surface<OrdinalType,ScalarType> sphere( 0,
						  1, 1, 1,
						  20, 20, 20,
						  10*10+10*10+10*10 - 10*10 );
  
  TEST_ASSERT( !sphere.isPlanar() );

  // Planar surface
  FACEMC::Surface<OrdinalType,ScalarType> plane( 0,
						 1, 2, 3,
						 4 );
  
  TEST_ASSERT( plane.isPlanar() );
}

UNIT_TEST_INSTANTIATION( Surface, isPlanar );

//---------------------------------------------------------------------------//
// Check that a surface can return the correct sense of a point
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Surface, 
				   getSense,
				   OrdinalType,
				   ScalarType )
{
  // Spherical surface centered at (-10,-10,-10) with radius 10
  FACEMC::Surface<OrdinalType,ScalarType> sphere( 0,
						  1, 1, 1,
						  20, 20, 20,
						  10*10+10*10+10*10 - 10*10 );
  
  FACEMC::Vector<ScalarType> point_1( 0.0, -10.0, -10.0 );
  FACEMC::Vector<ScalarType> point_2( -10.0, -10.0, -10.0 ); 
  FACEMC::Vector<ScalarType> point_3( 10.0, 10.0, 10.0 ); 

  FACEMC::SurfaceSense sense_1 = sphere.getSenseOfPoint( point_1 );
  FACEMC::SurfaceSense sense_2 = sphere.getSenseOfPoint( point_2 );
  FACEMC::SurfaceSense sense_3 = sphere.getSenseOfPoint( point_3 );

  TEST_EQUALITY_CONST( sense_1, FACEMC::ON_SURFACE );
  TEST_EQUALITY_CONST( sense_2, FACEMC::NEG_SURFACE_SENSE );
  TEST_EQUALITY_CONST( sense_3, FACEMC::POS_SURFACE_SENSE );
}

UNIT_TEST_INSTANTIATION( Surface, getSense );

//---------------------------------------------------------------------------//
// Check that a surface can return the correct unit normal at a point
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Surface, 
				   getUnitNormal,
				   OrdinalType,
				   ScalarType )
{
  // Spherical surface centered at (-10,-10,-10) with radius 10
  FACEMC::Surface<OrdinalType,ScalarType> sphere( 0,
						  1, 1, 1,
						  20, 20, 20,
						  10*10+10*10+10*10 - 10*10 );

  FACEMC::Vector<ScalarType> point_1( 0.0, -10.0, -10.0 ); 
  FACEMC::Vector<ScalarType> point_2( -10.0, 0.0, -10.0 );
  FACEMC::Vector<ScalarType> point_3( -10.0, -10.0, 0.0 );
  FACEMC::Vector<ScalarType> point_4( -20.0, -10.0, -10.0 );
  FACEMC::Vector<ScalarType> point_5( -10.0, -20.0, -10.0 );
  FACEMC::Vector<ScalarType> point_6( -10.0, -10.0, -20.0 );

  FACEMC::Vector<ScalarType> normal_1 = sphere.getUnitNormalAtPoint( point_1 );
  FACEMC::Vector<ScalarType> normal_2 = sphere.getUnitNormalAtPoint( point_2 );
  FACEMC::Vector<ScalarType> normal_3 = sphere.getUnitNormalAtPoint( point_3 );
  FACEMC::Vector<ScalarType> normal_4 = sphere.getUnitNormalAtPoint( point_4 );
  FACEMC::Vector<ScalarType> normal_5 = sphere.getUnitNormalAtPoint( point_5 );
  FACEMC::Vector<ScalarType> normal_6 = sphere.getUnitNormalAtPoint( point_6 );

  FACEMC::Vector<ScalarType> ref_normal_1( 1.0, 0.0, 0.0 );
  FACEMC::Vector<ScalarType> ref_normal_2( 0.0, 1.0, 0.0 );
  FACEMC::Vector<ScalarType> ref_normal_3( 0.0, 0.0, 1.0 );
  FACEMC::Vector<ScalarType> ref_normal_4( -1.0, 0.0, 0.0 );
  FACEMC::Vector<ScalarType> ref_normal_5( 0.0, -1.0, 0.0 );
  FACEMC::Vector<ScalarType> ref_normal_6( 0.0, 0.0, -1.0 );

  TEST_EQUALITY( normal_1, ref_normal_1 );
  TEST_EQUALITY( normal_2, ref_normal_2 );
  TEST_EQUALITY( normal_3, ref_normal_3 );
  TEST_EQUALITY( normal_4, ref_normal_4 );
  TEST_EQUALITY( normal_5, ref_normal_5 );
  TEST_EQUALITY( normal_6, ref_normal_6 );
}

UNIT_TEST_INSTANTIATION( Surface, getUnitNormal );

//---------------------------------------------------------------------------//
// Check that a surface can be created by transforming another surface with
// a translation vector.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Surface,
				   constructor_with_translation_vector,
				   OrdinalType,
				   ScalarType )
{
  FACEMC::Surface<OrdinalType,ScalarType> xy_plane( 0,
						    0, 0, 1,
						    0 );
  FACEMC::Surface<OrdinalType,ScalarType> unit_sphere( 1,
						       1, 1, 1,
						       0, 0, 0,
						       -1 );

    
  FACEMC::Vector<ScalarType> xaxis_translation_vector( 1.0, 0.0, 0.0 );

  FACEMC::Vector<ScalarType> yaxis_translation_vector( 0.0, 1.0, 0.0 );
  
  FACEMC::Vector<ScalarType> zaxis_translation_vector( 0.0, 0.0, 1.0 );

  FACEMC::Matrix<ScalarType> zero_matrix;

  FACEMC::Matrix<ScalarType> identity_matrix;
  identity_matrix.identity();
  
  // Create a plane at z=1 from the xy_plane
  FACEMC::Surface<OrdinalType,ScalarType> z_plane( 2,
						   xy_plane,
						   zaxis_translation_vector );

  FACEMC::Matrix<ScalarType> quad_form_matrix = 
    z_plane.getQuadraticFormMatrix();
  FACEMC::Vector<ScalarType> linear_term_vector = 
    z_plane.getLinearTermVector();

  FACEMC::Vector<ScalarType> ref_linear_term_vector( 0.0, 0.0, 1.0 );
  FACEMC::Matrix<ScalarType> ref_quad_form_matrix;
  ScalarType ref_constant_term = 1.0;

  TEST_EQUALITY( quad_form_matrix, ref_quad_form_matrix );
  TEST_COMPARE_FLOATING_ARRAYS( linear_term_vector(),
				ref_linear_term_vector(),
				Teuchos::ScalarTraits<ScalarType>::prec() );
  TEST_FLOATING_EQUALITY( z_plane.getConstantTerm(), 
			  ref_constant_term, 
			  Teuchos::ScalarTraits<ScalarType>::prec() );

  // Create a unit sphere centered at (-1, 0, 0) from the unit sphere at origin
  FACEMC::Surface<OrdinalType,ScalarType> unit_sphere_xaxis(
						    3,
						    unit_sphere,
						    xaxis_translation_vector );
  
  quad_form_matrix = unit_sphere_xaxis.getQuadraticFormMatrix();
  linear_term_vector = unit_sphere_xaxis.getLinearTermVector();

  ref_linear_term_vector =  FACEMC::createVector<ScalarType>( 2.0, 0.0, 0.0 );
  ref_quad_form_matrix.identity();
  ref_constant_term = 0.0;

  TEST_EQUALITY( quad_form_matrix, ref_quad_form_matrix );
  TEST_COMPARE_FLOATING_ARRAYS( linear_term_vector(),
				ref_linear_term_vector(),
				Teuchos::ScalarTraits<ScalarType>::prec() );
  TEST_FLOATING_EQUALITY( unit_sphere_xaxis.getConstantTerm(), 
			  ref_constant_term, 
			  Teuchos::ScalarTraits<ScalarType>::prec() );

  // Create a unit sphere centered at (0, -1, 0) from the unit sphere at origin
  FACEMC::Surface<OrdinalType,ScalarType> unit_sphere_yaxis(
						    3,
						    unit_sphere,
						    yaxis_translation_vector );

  quad_form_matrix = unit_sphere_yaxis.getQuadraticFormMatrix();
  linear_term_vector = unit_sphere_yaxis.getLinearTermVector();

  ref_linear_term_vector = FACEMC::createVector<ScalarType>( 0.0, 2.0, 0.0 );
  ref_quad_form_matrix.identity();
  ref_constant_term = 0.0;

  TEST_EQUALITY( quad_form_matrix, ref_quad_form_matrix );
  TEST_COMPARE_FLOATING_ARRAYS( linear_term_vector(),
				ref_linear_term_vector(),
				Teuchos::ScalarTraits<ScalarType>::prec() );
  TEST_FLOATING_EQUALITY( unit_sphere_yaxis.getConstantTerm(), 
			  ref_constant_term, 
			  Teuchos::ScalarTraits<ScalarType>::prec() );

  // Create a unit sphere centered at (0, 0, -1) from the unit sphere at origin
  FACEMC::Surface<OrdinalType,ScalarType> unit_sphere_zaxis( 
						    3,
						    unit_sphere,
						    zaxis_translation_vector );
  
  quad_form_matrix = unit_sphere_zaxis.getQuadraticFormMatrix();
  linear_term_vector = unit_sphere_zaxis.getLinearTermVector();

  ref_linear_term_vector = FACEMC::createVector<ScalarType>( 0.0, 0.0, 2.0 );
  ref_quad_form_matrix.identity();
  ref_constant_term = 0.0;

  TEST_EQUALITY( quad_form_matrix, ref_quad_form_matrix );
  TEST_COMPARE_FLOATING_ARRAYS( linear_term_vector(),
				ref_linear_term_vector(),
				Teuchos::ScalarTraits<ScalarType>::prec() );
  TEST_FLOATING_EQUALITY( unit_sphere_zaxis.getConstantTerm(), 
			  ref_constant_term, 
			  Teuchos::ScalarTraits<ScalarType>::prec() );
}

UNIT_TEST_INSTANTIATION( Surface, constructor_with_translation_vector );
  
//---------------------------------------------------------------------------//
// Check that a surface can be created by transforming another surface with
// a translation vector.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Surface,
				   constructor_with_rotation_matrix,
				   OrdinalType,
				   ScalarType )
{
  // General plane.
  FACEMC::Surface<OrdinalType,ScalarType> plane( 0,
						 5, 4, 3,
						 0 );
  // Unit skew cylinder with axis direction (1/2,1/2,sqrt(2)/2)
  FACEMC::Surface<OrdinalType,ScalarType> skew_cylinder( 1,
							 0.5, 1.0, 0.5,
							 0.0, 0.0, 1.0,
							 0.0, 0.0, 0.0,
							 -1.0 );

  // create the xy-plane by rotating the general plane
  FACEMC::Vector<ScalarType> desired_direction( 0.0, 0.0, 1.0 );
  FACEMC::Vector<ScalarType> current_direction( 1.0/2, 1.0/2, 1.0/sqrt(2.0) );
  FACEMC::Matrix<ScalarType> rotation_matrix =
    FACEMC::createRotationMatrixFromUnitVectors( desired_direction,
						 current_direction );

  FACEMC::Surface<OrdinalType,ScalarType> xy_plane( 2,
						    plane,
						    rotation_matrix );
  
  FACEMC::Matrix<ScalarType> quad_form_matrix = 
    xy_plane.getQuadraticFormMatrix();
  FACEMC::Vector<ScalarType> linear_term_vector = 
    xy_plane.getLinearTermVector();
  
  FACEMC::Matrix<ScalarType> ref_quad_form_matrix;
  FACEMC::Vector<ScalarType> ref_linear_term_vector( 0.0, 0.0, 1.0 );
  ScalarType ref_constant_term = 0.0;

  TEST_EQUALITY( quad_form_matrix, ref_quad_form_matrix );
  FACEMC_TEST_COMPARE_FLOATING_ARRAYS(
				    linear_term_vector(),
				    ref_linear_term_vector(),
				    Teuchos::ScalarTraits<ScalarType>::prec());
  TEST_FLOATING_EQUALITY( xy_plane.getConstantTerm(), 
			  ref_constant_term, 
			  Teuchos::ScalarTraits<ScalarType>::prec() );

  // create a cylinder parallel to z-axis centered at (1.0,1.0)
  FACEMC::Vector<ScalarType> eigenvalues =
  FACEMC::LinearAlgebra::computeEigenvaluesAndEigenvectors( 
					skew_cylinder.getQuadraticFormMatrix(),
  					rotation_matrix );
  if( !eigenvalues.hasUniqueElements() )
    FACEMC::LinearAlgebra::realignEigenvectors( eigenvalues, rotation_matrix );

  FACEMC::Surface<OrdinalType,ScalarType> z_cylinder( 3,
						      skew_cylinder,
						      rotation_matrix );

  quad_form_matrix = z_cylinder.getQuadraticFormMatrix();
  linear_term_vector = z_cylinder.getLinearTermVector();

  ref_quad_form_matrix = FACEMC::createMatrix<ScalarType>( 1.0,
							   0.0, 1.0,
							   0.0, 0.0, 0.0 );
  
  ref_linear_term_vector = FACEMC::createVector<ScalarType>( 0.0, 0.0, 0.0 );
  ref_constant_term = -1.0;

  FACEMC_TEST_COMPARE_FLOATING_ARRAYS(
				    quad_form_matrix(),
				    ref_quad_form_matrix(),
				    Teuchos::ScalarTraits<ScalarType>::prec());
  FACEMC_TEST_COMPARE_FLOATING_ARRAYS(
				    linear_term_vector(),
				    ref_linear_term_vector(),
				    Teuchos::ScalarTraits<ScalarType>::prec());
  TEST_FLOATING_EQUALITY( z_cylinder.getConstantTerm(), 
			  ref_constant_term, 
			  Teuchos::ScalarTraits<ScalarType>::prec() );
}

UNIT_TEST_INSTANTIATION( Surface, constructor_with_rotation_matrix );

//---------------------------------------------------------------------------//
// Check that a surface can be transformed with a rotation matrix and 
// a translation vector
// TEUCHOS_UNIT_TEST( Surface, transformSurface )
// {
//   // Spherical surface centered at (-10,-10,-10) with radius 10
//   Teuchos::RCP<FACEMC::Surface> surface_ptr( 
// 			     new FACEMC::Surface( 0,
// 						  1, 1, 1,
// 						  20, 20, 20,
// 						  10*10+10*10+10*10 - 10*10 ) 
// 					  );
  
//   // Direction vector of sphere local z-axis (in global coordinate system)
//   FACEMC::Surface::Vector initial_direction =
//     Teuchos::tuple( -sqrt(3.0)/3, -sqrt(3.0)/3, -sqrt(3.0)/3 );

//   // Direction vector of sphere local z-axis after transform
//   FACEMC::Surface::Vector final_direction = 
//     Teuchos::tuple( 0.0, 0.0, 1.0 );

//   // Rotation matrix
//   FACEMC::Surface::Matrix rotation_matrix =
//     FACEMC::LinearAlgebra::generateRotationMatrixFromUnitVectors(
// 				              initial_direction,
// 				              final_direction );

//   // Translation vector (vector from origin to sphere center)
//   FACEMC::Surface::Vector translation_vector = 
//     Teuchos::tuple( -10.0, -10.0, -10.0 );

//   surface_ptr->transformSurface( rotation_matrix,
// 				 translation_vector );

//   FACEMC::Surface::Matrix ref_quadratic_form_matrix = 
//     Teuchos::tuple( 1.0, 0.0, 0.0,
// 		    0.0, 1.0, 0.0,
// 		    0.0, 0.0, 1.0 );
//   FACEMC::Surface::Vector ref_linear_term_vector = 
//     Teuchos::tuple( 0.0, 0.0, 0.0 );
//   double ref_constant_term = -10*10;

//   FACEMC_TEST_COMPARE_FLOATING_ARRAYS( surface_ptr->getQuadraticFormMatrix()(), 
// 				       ref_quadratic_form_matrix(),
// 				       TOL );
//   FACEMC_TEST_COMPARE_FLOATING_ARRAYS( surface_ptr->getLinearTermVector()(), 
// 				       ref_linear_term_vector(),
// 				       TOL );
//   TEST_FLOATING_EQUALITY( surface_ptr->getConstantTerm(), 
// 			  ref_constant_term,
// 			  TOL );

//   // Spherical surface centered at (-10,-10,-10) with radius 10
//   surface_ptr.reset( new FACEMC::Surface( 0,
// 					  1, 1, 1,
// 					  20, 20, 20,
// 					  10*10+10*10+10*10 - 10*10 ) 
// 		     );
  
//   // Direction vector of sphere local z-axis (in global coordinate system)
//   initial_direction =
//     Teuchos::tuple( -sqrt(3.0)/3, -sqrt(3.0)/3, -sqrt(3.0)/3 );

//   // Direction vector of sphere local z-axis after transform
//   final_direction = 
//     Teuchos::tuple( 0.0, 1.0, 0.0 );

//   // Rotation matrix
//   rotation_matrix =
//     FACEMC::LinearAlgebra::generateRotationMatrixFromUnitVectors(
// 				              initial_direction,
// 					      final_direction );

//   // Translation vector 
//   translation_vector = 
//     Teuchos::tuple( -10.0, -10.0, -10.0 );

//   surface_ptr->transformSurface( rotation_matrix,
// 				 translation_vector );

//   translation_vector = 
//     Teuchos::tuple( 10.0, -10.0, -10.0 );

//   surface_ptr->transformSurface( ref_quadratic_form_matrix,
// 				 translation_vector );

//   ref_quadratic_form_matrix = 
//     Teuchos::tuple( 1.0, 0.0, 0.0,
// 		    0.0, 1.0, 0.0,
// 		    0.0, 0.0, 1.0 );
//   ref_linear_term_vector = 
//     Teuchos::tuple( 20.0, -20.0, -20.0 );
  
//   ref_constant_term = 10*10+10*10+10*10 - 10*10;

//   FACEMC_TEST_COMPARE_FLOATING_ARRAYS( surface_ptr->getQuadraticFormMatrix()(), 
// 				       ref_quadratic_form_matrix(),
// 				       TOL );
//   FACEMC_TEST_COMPARE_FLOATING_ARRAYS( surface_ptr->getLinearTermVector()(), 
// 				       ref_linear_term_vector(),
// 				       TOL );
//   TEST_FLOATING_EQUALITY( surface_ptr->getConstantTerm(), 
// 			  ref_constant_term,
// 			  TOL );

//   // cylindrical surface centered at (x0=-10, y0=-10) with radius 10
//   surface_ptr.reset( 
// 		    new FACEMC::Surface( 0,
// 					 1, 1, 0,
// 					 20, 20, 0,
// 					 10*10+10*10 - 10*10 ) 
// 		     );
  
//   // Direction vector of cylinder local z-axis (in global coordinate system)
//   initial_direction =
//     Teuchos::tuple( 0.0, 0.0, 1.0 );

//   // Direction vector of cylinder local z-axis after transformation
//   final_direction =
//     Teuchos::tuple( 1.0, 0.0, 0.0 );

//   // Rotation matrix
//   rotation_matrix =
//     FACEMC::LinearAlgebra::generateRotationMatrixFromUnitVectors(
// 				              initial_direction,
// 				              final_direction );

//   // Translation vector (vector from origin to cylinder center)
//   translation_vector = 
//     Teuchos::tuple( -10.0, -10.0, 0.0 );

//   surface_ptr->transformSurface( rotation_matrix,
// 				 translation_vector );

//   ref_quadratic_form_matrix = 
//     Teuchos::tuple( 0.0, 0.0, 0.0,
// 		    0.0, 1.0, 0.0,
// 		    0.0, 0.0, 1.0 );
//   ref_linear_term_vector = 
//     Teuchos::tuple( 0.0, 0.0, 0.0 );
//   ref_constant_term = -10*10;

//   FACEMC_TEST_COMPARE_FLOATING_ARRAYS( surface_ptr->getQuadraticFormMatrix()(), 
// 				       ref_quadratic_form_matrix(),
// 				       TOL );
//   FACEMC_TEST_COMPARE_FLOATING_ARRAYS( surface_ptr->getLinearTermVector()(), 
// 				       ref_linear_term_vector(),
// 				       TOL );
//   TEST_FLOATING_EQUALITY( surface_ptr->getConstantTerm(), 
// 			  ref_constant_term,
// 			  TOL );

//   // cylindrical surface centered at (x0=0, y0=-10) with radius 10
//   surface_ptr.reset( 
// 		    new FACEMC::Surface( 0,
// 					 1, 1, 0,
// 					 0, 20, 0,
// 					 10*10 - 10*10 ) 
// 		     );
  
//   // Direction vector of cylinder local z-axis (in global coordinate system)
//   initial_direction =
//     Teuchos::tuple( 0.0, 0.0, 1.0 );

//   // Direction vector of cylinder local z-axis after transformation
//   final_direction =
//     Teuchos::tuple( 1.0, 0.0, 0.0 );

//   // Rotation matrix
//   rotation_matrix =
//     FACEMC::LinearAlgebra::generateRotationMatrixFromUnitVectors(
// 				              initial_direction,
// 				              final_direction );

//   // Translation vector (vector from origin to cylinder center)
//   translation_vector = 
//     Teuchos::tuple( 0.0, -10.0, 0.0 );

//   surface_ptr->transformSurface( rotation_matrix,
// 				 translation_vector );

//   ref_quadratic_form_matrix = 
//     Teuchos::tuple( 0.0, 0.0, 0.0,
// 		    0.0, 1.0, 0.0,
// 		    0.0, 0.0, 1.0 );
//   ref_linear_term_vector = 
//     Teuchos::tuple( 0.0, 0.0, 0.0 );
//   ref_constant_term = -10*10;

//   FACEMC_TEST_COMPARE_FLOATING_ARRAYS( surface_ptr->getQuadraticFormMatrix()(), 
// 				       ref_quadratic_form_matrix(),
// 				       TOL );
//   FACEMC_TEST_COMPARE_FLOATING_ARRAYS( surface_ptr->getLinearTermVector()(), 
// 				       ref_linear_term_vector(),
// 				       TOL );
//   TEST_FLOATING_EQUALITY( surface_ptr->getConstantTerm(), 
// 			  ref_constant_term,
// 			  TOL );

//   // planar surface 
//   surface_ptr.reset( 
// 		    new FACEMC::Surface( 0,
// 					 6, 5, 4,
// 					 3 ) 
// 		     );
  
//   // Direction vector of plane local z-axis (in global coordinate system)
//   initial_direction =
//     Teuchos::tuple( 6.0/sqrt(6.0*6.0+5.0*5.0+4.0*4.0), 
// 		    5.0/sqrt(6.0*6.0+5.0*5.0+4.0*4.0), 
// 		    4.0/sqrt(6.0*6.0+5.0*5.0+4.0*4.0) );

//   // Direction vector of plane local z-axis after transformation
//   final_direction =
//     Teuchos::tuple( 0.0, 0.0, 1.0 );

//   // Rotation matrix
//   rotation_matrix =
//     FACEMC::LinearAlgebra::generateRotationMatrixFromUnitVectors(
// 					      initial_direction,
// 				              final_direction );
  
//   rotation_matrix = FACEMC::LinearAlgebra::matrixTranspose( rotation_matrix );
  
//   // Translation vector (vector from origin to plane)
//   translation_vector = 
//     Teuchos::tuple( 0.0, 0.0, surface_ptr->getConstantTerm()/-4.0 );

//   surface_ptr->transformSurface( rotation_matrix,
// 				 translation_vector );

//   ref_quadratic_form_matrix = 
//     Teuchos::tuple( 0.0, 0.0, 0.0,
// 		    0.0, 0.0, 0.0,
// 		    0.0, 0.0, 0.0 );
//   ref_linear_term_vector = 
//     Teuchos::tuple( 0.0, 0.0, sqrt(6.0*6.0+5.0*5.0+4.0*4.0) );
//   ref_constant_term = 0.0;

//   FACEMC_TEST_COMPARE_FLOATING_ARRAYS( surface_ptr->getQuadraticFormMatrix()(), 
// 				       ref_quadratic_form_matrix(),
// 				       TOL );
//   FACEMC_TEST_COMPARE_FLOATING_ARRAYS( surface_ptr->getLinearTermVector()(), 
// 				       ref_linear_term_vector(),
// 				       TOL );
//   TEST_FLOATING_EQUALITY( surface_ptr->getConstantTerm(), 
// 			  ref_constant_term,
// 			  TOL );

// }  
	   
//---------------------------------------------------------------------------//
// end tstSurface.cpp
//---------------------------------------------------------------------------//
