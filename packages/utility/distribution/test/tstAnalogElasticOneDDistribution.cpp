//---------------------------------------------------------------------------//
//!
//! \file   tstAnalogElasticOneDDistribution.cpp
//! \author Luke Kersting
//! \brief  Analog elastic distribution unit tests.
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Boost Includes
#include <boost/units/systems/si.hpp>
#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/cgs.hpp>
#include <boost/units/io.hpp>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_VerboseObject.hpp>

// FRENSIE Includes
#include "Utility_UnitTestHarnessExtensions.hpp"
#include "Utility_TabularOneDDistribution.hpp"
#include "Utility_AnalogElasticOneDDistribution.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_PhysicalConstants.hpp"
#include "Utility_UnitTraits.hpp"
#include "Utility_QuantityTraits.hpp"

using boost::units::quantity;
using namespace Utility::Units;
namespace si = boost::units::si;
namespace cgs = boost::units::cgs;

typedef quantity<si::dimensionless> dl;

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//

Teuchos::RCP<Teuchos::ParameterList> test_dists_list;

Teuchos::RCP<Utility::OneDDistribution> distribution;
Teuchos::RCP<Utility::TabularOneDDistribution> tab_distribution;

Teuchos::RCP<Utility::UnitAwareOneDDistribution<si::dimensionless,si::amount> >
  unit_aware_distribution;
Teuchos::RCP<Utility::UnitAwareTabularOneDDistribution<si::dimensionless,si::amount> >
unit_aware_tab_distribution;

double eta = 1e-12;
double cross_section_ratio = 0.1;

//---------------------------------------------------------------------------//
// Instantiation Macros.
//---------------------------------------------------------------------------//
#define UNIT_TEST_INSTANTIATION( type, name )                                \
  typedef Utility::LinLin LinLin;                                        \
  typedef Utility::LogLin LogLin;                                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( type, name, LinLin )                \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( type, name, LogLin )

//---------------------------------------------------------------------------//
// Testing Functions.
//---------------------------------------------------------------------------//
// Initialize the distribution
template<typename InterpolationPolicy, typename BaseDistribution>
void initialize( Teuchos::RCP<BaseDistribution>& dist )
{
  // Use the basic constructor
  Teuchos::Array<typename BaseDistribution::IndepQuantity>
    independent_values( 4 );
  Utility::setQuantity( independent_values[0], -1.0 );
  Utility::setQuantity( independent_values[1], 0.0 );
  Utility::setQuantity( independent_values[2], 0.5 );
  Utility::setQuantity( independent_values[3], 0.999999 );

  Teuchos::Array<typename BaseDistribution::DepQuantity> dependent_values( 4 );
  Utility::setQuantity( dependent_values[0], 1e2 );
  Utility::setQuantity( dependent_values[1], 1e1 );
  Utility::setQuantity( dependent_values[2], 1.0 );
  Utility::setQuantity( dependent_values[3], 1e-1 );

  dist.reset(new Utility::UnitAwareAnalogElasticOneDDistribution<InterpolationPolicy,typename BaseDistribution::IndepUnit, typename BaseDistribution::DepUnit>(
                                                      independent_values,
                                                      dependent_values,
                                                      eta,
                                                      cross_section_ratio ) );
}

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that the distribution can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   evaluate,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_EQUALITY_CONST( distribution->evaluate( -1.0 ), 1e1 );
  TEST_EQUALITY_CONST( distribution->evaluate( 0.0 ), 1.0 );
  TEST_EQUALITY_CONST( distribution->evaluate( 0.5 ), 1e-1 );
  TEST_FLOATING_EQUALITY( distribution->evaluate( 0.999999 ), 1e-2, 1e-15 );
  TEST_FLOATING_EQUALITY( distribution->evaluate( 1.0 ), 1.00000200000100021e+10, 1e-15 );

  double sample = distribution->evaluate( distribution->sample() );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, evaluate );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   evaluate,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_EQUALITY_CONST( unit_aware_distribution->evaluate( -1.0 ),
                       1e1*si::mole );
  TEST_EQUALITY_CONST( unit_aware_distribution->evaluate( 0.0 ),
                       1.0*si::mole );
  TEST_EQUALITY_CONST( unit_aware_distribution->evaluate( 0.5 ),
                       1e-1*si::mole );
  UTILITY_TEST_FLOATING_EQUALITY( unit_aware_distribution->evaluate( 0.999999 ),
                                  1e-2*si::mole,
                                  1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY( unit_aware_distribution->evaluate( 1.0 ),
                                  1.00000200000100021e+10*si::mole,
                                  1e-15 );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution, evaluate );

//---------------------------------------------------------------------------//
// Check that the PDF can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   evaluatePDF,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( -1.0 ),
                          1.7233951046594968293e-1,
                          1e-6 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 0.0 ),
                          1.7233951046594968293e-2,
                          1e-6 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 0.5 ),
                          1.7233951046594968293e-3,
                          1e-6 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 0.999999 ),
                          1.7233951046594968293e-4,
                          1e-6 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 1.0 ),
                          1.7233985514514300972e9,
                          1e-6 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, evaluatePDF );

//---------------------------------------------------------------------------//
// Check that the unit-aware PDF can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   evaluatePDF,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  UTILITY_TEST_FLOATING_EQUALITY(
                             unit_aware_distribution->evaluatePDF( -1.0*si::dimensionless() ),
                             1.7233951046594968293e-1/si::dimensionless(),
                             1e-6 );
  UTILITY_TEST_FLOATING_EQUALITY(
                             unit_aware_distribution->evaluatePDF( 0.0*si::dimensionless() ),
                             1.7233951046594968293e-2/si::dimensionless(),
                             1e-6 );
  UTILITY_TEST_FLOATING_EQUALITY(
                             unit_aware_distribution->evaluatePDF( 0.5*si::dimensionless() ),
                             1.7233951046594968293e-3/si::dimensionless(),
                             1e-6 );
  UTILITY_TEST_FLOATING_EQUALITY(
                              unit_aware_distribution->evaluatePDF( 0.999999*si::dimensionless() ),
                              1.7233951046594968293e-4/si::dimensionless(),
                              1e-6 );
  UTILITY_TEST_FLOATING_EQUALITY(
                              unit_aware_distribution->evaluatePDF( 1.0*si::dimensionless() ),
                              1.7233985514514300972e9/si::dimensionless(),
                              1e-6 );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution, evaluatePDF );

//---------------------------------------------------------------------------//
// Check that the CDF can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   evaluateCDF,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( tab_distribution );

  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( -1.0 ),
                          0.0000000000,
                          1e-10 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 0.0 ),
                          0.094786730756272336018,
                          1e-10 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 0.5 ),
                          0.09952606729408594588,
                          1e-10 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 0.999999 ),
                          0.1,
                          1e-10 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 1.0 ),
                          1.0000000000,
                          1e-10 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, evaluateCDF );

//---------------------------------------------------------------------------//
// Check that the unit-aware CDF can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   evaluateCDF,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_tab_distribution );

  TEST_FLOATING_EQUALITY( unit_aware_tab_distribution->evaluateCDF( -1.0*si::dimensionless() ),
                          0.0000000000,
                          1e-10 );
  TEST_FLOATING_EQUALITY( unit_aware_tab_distribution->evaluateCDF( 0.0*si::dimensionless() ),
                          0.094786730756272336018,
                          1e-10 );
  TEST_FLOATING_EQUALITY( unit_aware_tab_distribution->evaluateCDF( 0.5*si::dimensionless() ),
                          0.09952606729408594588,
                          1e-10 );
  TEST_FLOATING_EQUALITY( unit_aware_tab_distribution->evaluateCDF( 0.999999*si::dimensionless() ),
                          0.1,
                          1e-10 );
  TEST_FLOATING_EQUALITY( unit_aware_tab_distribution->evaluateCDF( 1.0*si::dimensionless() ),
                          1.0,
                          1e-10 );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution, evaluateCDF );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   sample,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  std::vector<double> fake_stream( 3 );
  fake_stream[0] = 0.0;
  fake_stream[1] = cross_section_ratio;
  fake_stream[2] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double sample = distribution->sample();
  TEST_EQUALITY_CONST( sample, -1.0 );

  sample = distribution->sample();
  TEST_FLOATING_EQUALITY( sample, 0.999999, 1e-12 );

  sample = distribution->sample();
  TEST_FLOATING_EQUALITY( sample, 1.0, 1e-12 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = distribution->sample();
  TEST_COMPARE( sample, >=, -1.0 );
  TEST_COMPARE( sample, <=, 1.0 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, sample );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   sample,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  std::vector<double> fake_stream( 3 );
  fake_stream[0] = 0.0;
  fake_stream[1] = cross_section_ratio;
  fake_stream[2] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  quantity<si::dimensionless> sample = unit_aware_distribution->sample();
  TEST_EQUALITY_CONST( sample, -1.0*si::dimensionless() );

  sample = unit_aware_distribution->sample();
  UTILITY_TEST_FLOATING_EQUALITY( sample, 0.999999*si::dimensionless(), 1e-12 );

  sample = unit_aware_distribution->sample();
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1.0*si::dimensionless(), 1e-12 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = unit_aware_distribution->sample();
  TEST_COMPARE( sample, >=, -1.0*si::dimensionless() );
  TEST_COMPARE( sample, <=, 1.0*si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution, sample );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   sampleAndRecordTrials,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  std::vector<double> fake_stream( 3 );
  fake_stream[0] = 0.0;
  fake_stream[1] = cross_section_ratio;
  fake_stream[2] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  unsigned trials = 0;

  double sample = distribution->sampleAndRecordTrials( trials );
  TEST_EQUALITY_CONST( sample, -1.0 );
  TEST_EQUALITY_CONST( 1.0/trials, 1.0 );

  sample = distribution->sampleAndRecordTrials( trials );
  TEST_FLOATING_EQUALITY( sample, 0.999999, 1e-12 );
  TEST_EQUALITY_CONST( 2.0/trials, 1.0 );

  sample = distribution->sampleAndRecordTrials( trials );
  TEST_FLOATING_EQUALITY( sample, 1.0, 1e-12 );
  TEST_EQUALITY_CONST( 3.0/trials, 1.0 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = distribution->sampleAndRecordTrials( trials );
  TEST_COMPARE( sample, >=, -1.0 );
  TEST_COMPARE( sample, <=, 1.0 );
  TEST_EQUALITY_CONST( 4.0/trials, 1.0 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, sampleAndRecordTrials );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   sampleAndRecordTrials,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  std::vector<double> fake_stream( 3 );
  fake_stream[0] = 0.0;
  fake_stream[1] = cross_section_ratio;
  fake_stream[2] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  unsigned trials = 0;

  quantity<si::dimensionless> sample =
    unit_aware_distribution->sampleAndRecordTrials( trials );
  TEST_EQUALITY_CONST( sample, -1.0*si::dimensionless() );
  TEST_EQUALITY_CONST( 1.0/trials, 1.0 );

  sample = unit_aware_distribution->sampleAndRecordTrials( trials );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 0.999999*si::dimensionless(), 1e-12 );
  TEST_EQUALITY_CONST( 2.0/trials, 1.0 );

  sample = unit_aware_distribution->sampleAndRecordTrials( trials );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1.0*si::dimensionless(), 1e-12 );
  TEST_EQUALITY_CONST( 3.0/trials, 1.0 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = unit_aware_distribution->sampleAndRecordTrials( trials );
  TEST_COMPARE( sample, >=, -1.0*si::dimensionless() );
  TEST_COMPARE( sample, <=, 1.0*si::dimensionless() );
  TEST_EQUALITY_CONST( 4.0/trials, 1.0 );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution, sampleAndRecordTrials );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   sampleAndRecordBinIndex,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( tab_distribution );

  std::vector<double> fake_stream( 3 );
  fake_stream[0] = 0.0;
  fake_stream[1] = cross_section_ratio;
  fake_stream[2] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  unsigned bin_index;

  double sample = tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_EQUALITY_CONST( sample, -1.0 );
  TEST_EQUALITY_CONST( bin_index, 0u );

  sample = tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_FLOATING_EQUALITY( sample, 0.999999, 1e-12 );
  TEST_EQUALITY_CONST( bin_index, 3u );

  sample = tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_FLOATING_EQUALITY( sample, 1.0, 1e-12 );
  TEST_EQUALITY_CONST( bin_index, 3u );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_COMPARE( sample, >=, -1.0 );
  TEST_COMPARE( sample, <=, 1.0 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, sampleAndRecordBinIndex );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   sampleAndRecordBinIndex,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_tab_distribution );

  std::vector<double> fake_stream( 3 );
  fake_stream[0] = 0.0;
  fake_stream[1] = cross_section_ratio;
  fake_stream[2] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  unsigned bin_index;

  quantity<si::dimensionless> sample =
    unit_aware_tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_EQUALITY_CONST( sample, -1.0*si::dimensionless() );
  TEST_EQUALITY_CONST( bin_index, 0u );

  sample = unit_aware_tab_distribution->sampleAndRecordBinIndex( bin_index );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 0.999999*si::dimensionless(), 1e-12 );
  TEST_EQUALITY_CONST( bin_index, 3u );

  sample = unit_aware_tab_distribution->sampleAndRecordBinIndex( bin_index );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1.0*si::dimensionless(), 1e-12 );
  TEST_EQUALITY_CONST( bin_index, 3u );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = unit_aware_tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_COMPARE( sample, >=, -1.0*si::dimensionless() );
  TEST_COMPARE( sample, <=, 1.0*si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution,
                         sampleAndRecordBinIndex );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   sampleWithRandomNumber,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( tab_distribution );

  double sample = tab_distribution->sampleWithRandomNumber( 0.0 );
  TEST_EQUALITY_CONST( sample, -1.0 );

  sample = tab_distribution->sampleWithRandomNumber( cross_section_ratio );
  TEST_FLOATING_EQUALITY( sample, 0.999999, 1e-12 );

  sample = tab_distribution->sampleWithRandomNumber( 1.0 - 1e-15 );
  TEST_FLOATING_EQUALITY( sample, 1.0, 1e-12 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, sampleWithRandomNumber );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   sampleWithRandomNumber,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_tab_distribution );

  quantity<si::dimensionless> sample =
    unit_aware_tab_distribution->sampleWithRandomNumber( 0.0 );
  TEST_EQUALITY_CONST( sample, -1.0*si::dimensionless() );

  sample = unit_aware_tab_distribution->sampleWithRandomNumber( cross_section_ratio );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 0.999999*si::dimensionless(), 1e-12 );

  sample = unit_aware_tab_distribution->sampleWithRandomNumber( 1.0 - 1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1.0*si::dimensionless(), 1e-12 );

  double random_number = cross_section_ratio/2.0;
  sample = unit_aware_tab_distribution->sampleWithRandomNumber( random_number );
  TEST_COMPARE( sample, >, -1.0 );
  TEST_COMPARE( sample, <, 0.999999 );

  random_number = (1.0 - cross_section_ratio)/2.0;
  sample = unit_aware_tab_distribution->sampleWithRandomNumber( random_number );
  TEST_COMPARE( sample, >, 0.999999 );
  TEST_COMPARE( sample, <, 1.0 );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution,
                         sampleWithRandomNumber );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled from a subrange
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   sampleInSubrange,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( tab_distribution );

  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 1.0 - 1e-15;
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double sample = tab_distribution->sampleInSubrange( 1e-1  );
  TEST_EQUALITY_CONST( sample, -1.0 );

  sample = tab_distribution->sampleInSubrange( 1e-1 );
  TEST_FLOATING_EQUALITY( sample, 1e-1, 1e-12 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = tab_distribution->sampleInSubrange( 1e-1 );
  TEST_COMPARE( sample, >=, -1.0 );
  TEST_COMPARE( sample, <=, 1e-1 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, sampleInSubrange );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled from a subrange
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   sampleInSubrange,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_tab_distribution );

  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  quantity<si::dimensionless> sample =
    unit_aware_tab_distribution->sampleInSubrange( 1e-1*si::dimensionless()  );
  TEST_EQUALITY_CONST( sample, -1.0*si::dimensionless() );

  sample = unit_aware_tab_distribution->sampleInSubrange( 1e-1*si::dimensionless() );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1e-1*si::dimensionless(), 1e-12 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = unit_aware_tab_distribution->sampleInSubrange( 1e-1*si::dimensionless() );
  TEST_COMPARE( sample, >=, -1.0*si::dimensionless() );
  TEST_COMPARE( sample, <=, 1e-1*si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution, sampleInSubrange );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled from a subrange
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   sampleWithRandomNumberInSubrange,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( tab_distribution );

  double sample =
    tab_distribution->sampleWithRandomNumberInSubrange( 0.0, 1e-1  );
  TEST_EQUALITY_CONST( sample, -1.0 );

  sample = tab_distribution->sampleWithRandomNumberInSubrange( 1.0, 1e-1 );
  TEST_FLOATING_EQUALITY( sample, 1e-1, 1e-12 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution,
                         sampleWithRandomNumberInSubrange );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled from a subrange
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   sampleWithRandomNumberInSubrange,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_tab_distribution );

  quantity<si::dimensionless> sample =
    unit_aware_tab_distribution->sampleWithRandomNumberInSubrange(
                                                               0.0, 1e-1*si::dimensionless() );
  TEST_EQUALITY_CONST( sample, -1.0*si::dimensionless() );

  sample = unit_aware_tab_distribution->sampleWithRandomNumberInSubrange(
                                                               1.0, 1e-1*si::dimensionless() );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1e-1*si::dimensionless(), 1e-12 );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution,
                         sampleWithRandomNumberInSubrange );

//---------------------------------------------------------------------------//
// Check that the upper bound of the distribution independent variable can be
// returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   getUpperBoundOfIndepVar,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_EQUALITY_CONST( distribution->getUpperBoundOfIndepVar(), 1.0 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, getUpperBoundOfIndepVar );

//---------------------------------------------------------------------------//
// Check that the upper bound of the unit-aware distribution independent
// variable can be returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   getUpperBoundOfIndepVar,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_EQUALITY_CONST( unit_aware_distribution->getUpperBoundOfIndepVar(),
                       1.0*si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution,
                         getUpperBoundOfIndepVar );

//---------------------------------------------------------------------------//
// Check that the cutoff bound of the distribution independent variable can be
// returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   getCutoffBoundOfIndepVar,
                                   InterpolationPolicy )
{
  Teuchos::RCP<Utility::AnalogElasticOneDDistribution<InterpolationPolicy> >
                analog_distribution;

  initialize<InterpolationPolicy>( analog_distribution );

  TEST_EQUALITY_CONST( analog_distribution->getCutoffBoundOfIndepVar(), 0.999999 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, getCutoffBoundOfIndepVar );

//---------------------------------------------------------------------------//
// Check that the cutoff bound of the unit-aware distribution independent
// variable can be returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   getCutoffBoundOfIndepVar,
                                   InterpolationPolicy )
{
Teuchos::RCP<Utility::UnitAwareAnalogElasticOneDDistribution<InterpolationPolicy,si::dimensionless,si::amount> >
unit_aware_analog_distribution;

  initialize<InterpolationPolicy>( unit_aware_analog_distribution );

  TEST_EQUALITY_CONST( unit_aware_analog_distribution->getCutoffBoundOfIndepVar(),
                       0.999999*si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution,
                         getCutoffBoundOfIndepVar );

//---------------------------------------------------------------------------//
// Check that the lower bound of the distribution independent variable can be
// returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   getLowerBoundOfIndepVar,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_EQUALITY_CONST( distribution->getLowerBoundOfIndepVar(), -1.0 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, getLowerBoundOfIndepVar );

//---------------------------------------------------------------------------//
// Check that the lower bound of the unit-aware distribution independent
// variable can be returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   getLowerBoundOfIndepVar,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_EQUALITY_CONST( unit_aware_distribution->getLowerBoundOfIndepVar(),
                       -1.0*si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution,
                         getLowerBoundOfIndepVar );

//---------------------------------------------------------------------------//
// Check that the distribution type can be returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   getDistributionType,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_EQUALITY_CONST( distribution->getDistributionType(),
                       Utility::ANALOG_ELASTIC_DISTRIBUTION );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, getDistributionType );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution type can be returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   getDistributionType,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_EQUALITY_CONST( unit_aware_distribution->getDistributionType(),
                       Utility::ANALOG_ELASTIC_DISTRIBUTION );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution, getDistributionType );

//---------------------------------------------------------------------------//
// Check if the distribution is tabular
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   isTabular,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_ASSERT( distribution->isTabular() );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, isTabular );

//---------------------------------------------------------------------------//
// Check if the unit-aware distribution is tabular
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   isTabular,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_ASSERT( unit_aware_distribution->isTabular() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution, isTabular );

//---------------------------------------------------------------------------//
// Check that the distribution is continuous
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   isContinuous,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_ASSERT( distribution->isContinuous() );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, isContinuous );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution is continuous
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   isContinuous,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_ASSERT( unit_aware_distribution->isContinuous() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution, isContinuous );

//---------------------------------------------------------------------------//
// Check if the distribution is compatible with the interpolation type
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   isCompatibleWithInterpType,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_ASSERT( distribution->isCompatibleWithInterpType<Utility::LinLin>() );
  TEST_ASSERT( distribution->isCompatibleWithInterpType<Utility::LinLog>() );
  TEST_ASSERT( distribution->isCompatibleWithInterpType<Utility::LogLin>() );
  TEST_ASSERT( distribution->isCompatibleWithInterpType<Utility::LogLog>() );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, isCompatibleWithInterpType );

//---------------------------------------------------------------------------//
// Check if the unit-aware distribution is compatible with the interp type
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   isCompatibleWithInterpType,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );


  TEST_ASSERT( unit_aware_distribution->isCompatibleWithInterpType<Utility::LinLin>() );
  TEST_ASSERT( unit_aware_distribution->isCompatibleWithInterpType<Utility::LinLog>() );
  TEST_ASSERT( unit_aware_distribution->isCompatibleWithInterpType<Utility::LogLin>() );
  TEST_ASSERT( unit_aware_distribution->isCompatibleWithInterpType<Utility::LogLog>() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution,
                         isCompatibleWithInterpType );

//---------------------------------------------------------------------------//
// Check that the distribution can be written to an xml file
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticOneDDistribution,
                                   toParameterList,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  typedef Utility::AnalogElasticOneDDistribution<InterpolationPolicy> Distribution;

  Teuchos::RCP<Distribution> true_distribution =
    Teuchos::rcp_dynamic_cast<Distribution>( distribution );

  Teuchos::ParameterList parameter_list;

  parameter_list.set<Distribution>( "test distribution",
                                      *true_distribution );

  std::ostringstream xml_file_name;
  xml_file_name << "analog_elastic_" << InterpolationPolicy::name()
                << "_dist_test_list.xml";

  Teuchos::writeParameterListToXmlFile( parameter_list,
                                        xml_file_name.str() );

  Teuchos::RCP<Teuchos::ParameterList> read_parameter_list =
    Teuchos::getParametersFromXmlFile( xml_file_name.str() );

  TEST_EQUALITY( parameter_list, *read_parameter_list );

  Teuchos::RCP<Distribution>
    copy_distribution( new Distribution );

  *copy_distribution =
        read_parameter_list->get<Distribution>( "test distribution");

  TEST_EQUALITY( *copy_distribution, *true_distribution );
}

UNIT_TEST_INSTANTIATION( AnalogElasticOneDDistribution, toParameterList );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be written to an xml file
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   toParameterList,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  typedef Utility::UnitAwareAnalogElasticOneDDistribution<InterpolationPolicy,si::dimensionless,si::amount> Distribution;

  Teuchos::RCP<Distribution> true_distribution =
    Teuchos::rcp_dynamic_cast<Distribution>( unit_aware_distribution );

  Teuchos::ParameterList parameter_list;

  parameter_list.set<Distribution>( "test distribution",
                                    *true_distribution );

  std::ostringstream xml_file_name;
  xml_file_name << "unit_aware_tabular_" << InterpolationPolicy::name()
                << "_dist_test_list.xml";

  Teuchos::writeParameterListToXmlFile( parameter_list,
                                        xml_file_name.str() );

  Teuchos::RCP<Teuchos::ParameterList> read_parameter_list =
    Teuchos::getParametersFromXmlFile( xml_file_name.str() );

  TEST_EQUALITY( parameter_list, *read_parameter_list );

  Teuchos::RCP<Distribution>
    copy_distribution( new Distribution );

  *copy_distribution = read_parameter_list->get<Distribution>(
                                                          "test distribution");

  TEST_EQUALITY( *copy_distribution, *true_distribution );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticOneDDistribution, toParameterList );

//---------------------------------------------------------------------------//
// Check that the distribution can be read from an xml file
TEUCHOS_UNIT_TEST( AnalogElasticOneDDistribution, fromParameterList )
{
  Utility::AnalogElasticOneDDistribution<Utility::LinLin> distribution_1 =
    test_dists_list->get<Utility::AnalogElasticOneDDistribution<Utility::LinLin> >( "Analog Elastic Distribution A" );

  TEST_EQUALITY_CONST( distribution_1.getLowerBoundOfIndepVar(), -1.0 );
  TEST_EQUALITY_CONST( distribution_1.getUpperBoundOfIndepVar(), 1.0 );
  TEST_EQUALITY_CONST( distribution_1.getMoliereScreeningConstant(), 1.0 );
  TEST_EQUALITY_CONST( distribution_1.getCutoffCrossSectionRatio(), 0.1 );

  distribution_1 =
    test_dists_list->get<Utility::AnalogElasticOneDDistribution<Utility::LinLin> >( "Analog Elastic Distribution B" );

  TEST_EQUALITY_CONST( distribution_1.getLowerBoundOfIndepVar(), -1.0 );
  TEST_EQUALITY_CONST( distribution_1.getUpperBoundOfIndepVar(), 1.0 );
  TEST_EQUALITY_CONST( distribution_1.getMoliereScreeningConstant(), 1.0 );
  TEST_EQUALITY_CONST( distribution_1.getCutoffCrossSectionRatio(), 0.2 );

  Utility::AnalogElasticOneDDistribution<Utility::LogLin> distribution_2 =
    test_dists_list->get<Utility::AnalogElasticOneDDistribution<Utility::LogLin> >( "Analog Elastic Distribution C" );

  TEST_EQUALITY_CONST( distribution_2.getLowerBoundOfIndepVar(), -1.0 );
  TEST_EQUALITY_CONST( distribution_2.getUpperBoundOfIndepVar(), 1.0 );
  TEST_EQUALITY_CONST( distribution_2.getMoliereScreeningConstant(), 1.0 );
  TEST_EQUALITY_CONST( distribution_2.getCutoffCrossSectionRatio(), 0.5 );
}

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be read from an xml file
TEUCHOS_UNIT_TEST( UnitAwareAnalogElasticOneDDistribution, fromParameterList )
{
  Utility::UnitAwareAnalogElasticOneDDistribution<Utility::LinLin,si::dimensionless,si::amount>
    distribution_1 =
    test_dists_list->get<Utility::UnitAwareAnalogElasticOneDDistribution<Utility::LinLin,si::dimensionless,si::amount> >( "Unit-Aware Analog Elastic Distribution A" );

  TEST_EQUALITY_CONST( distribution_1.getLowerBoundOfIndepVar(), -1.0 );
  TEST_EQUALITY_CONST( distribution_1.getUpperBoundOfIndepVar(), 1.0 );
  TEST_EQUALITY_CONST( distribution_1.getMoliereScreeningConstant(), 1.0 );
  TEST_EQUALITY_CONST( distribution_1.getCutoffCrossSectionRatio(), 0.1 );
  UTILITY_TEST_FLOATING_EQUALITY(
                        distribution_1.evaluateCDF( 0.99999911111120975971 ),
                        0.2,
                        1e-9 );
  TEST_EQUALITY_CONST( distribution_1.sampleWithRandomNumber( 0.2 ),
                       0.99999911111120975971 );

  distribution_1 =
    test_dists_list->get<Utility::UnitAwareAnalogElasticOneDDistribution<Utility::LinLin,si::dimensionless,si::amount> >( "Unit-Aware Analog Elastic Distribution B" );

  TEST_EQUALITY_CONST( distribution_1.getLowerBoundOfIndepVar(), -1.0 );
  TEST_EQUALITY_CONST( distribution_1.getUpperBoundOfIndepVar(), 1.0 );
  TEST_EQUALITY_CONST( distribution_1.getMoliereScreeningConstant(), 1.0 );
  TEST_EQUALITY_CONST( distribution_1.getCutoffCrossSectionRatio(), 0.2 );
  UTILITY_TEST_FLOATING_EQUALITY(
                        distribution_1.evaluateCDF( 0.99999912500010945671 ),
                        0.3,
                        1e-9 );
  TEST_EQUALITY_CONST( distribution_1.sampleWithRandomNumber( 0.3 ), 0.99999912500010945671 );

  Utility::UnitAwareAnalogElasticOneDDistribution<Utility::LogLin,si::dimensionless,si::amount>
    distribution_2 =
    test_dists_list->get<Utility::UnitAwareAnalogElasticOneDDistribution<Utility::LogLin,si::dimensionless,si::amount> >( "Unit-Aware Analog Elastic Distribution C" );

  TEST_EQUALITY_CONST( distribution_2.getLowerBoundOfIndepVar(), -1.0 );
  TEST_EQUALITY_CONST( distribution_2.getUpperBoundOfIndepVar(), 1.0 );
  TEST_EQUALITY_CONST( distribution_2.getMoliereScreeningConstant(), 1.0 );
  TEST_EQUALITY_CONST( distribution_2.getCutoffCrossSectionRatio(), 0.5 );
  TEST_EQUALITY_CONST( distribution_2.evaluateCDF( 0.9999995 ),
                       7.49999875035006980e-01 );
  TEST_EQUALITY_CONST( distribution_2.sampleWithRandomNumber( 7.49999875035006980e-01 ),
                       0.9999995 );
}

//---------------------------------------------------------------------------//
// Check that distributions can be scaled
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( UnitAwareAnalogElasticOneDDistribution,
                                   explicit_conversion,
                                   IndepUnitA,
                                   DepUnitA,
                                   IndepUnitB,
                                   DepUnitB )
{
  typedef typename Utility::UnitTraits<IndepUnitA>::template GetQuantityType<double>::type IndepQuantityA;
  typedef typename Utility::UnitTraits<typename Utility::UnitTraits<IndepUnitA>::InverseUnit>::template GetQuantityType<double>::type InverseIndepQuantityA;

  typedef typename Utility::UnitTraits<IndepUnitB>::template GetQuantityType<double>::type IndepQuantityB;
  typedef typename Utility::UnitTraits<typename Utility::UnitTraits<IndepUnitB>::InverseUnit>::template GetQuantityType<double>::type InverseIndepQuantityB;

  typedef typename Utility::UnitTraits<DepUnitA>::template GetQuantityType<double>::type DepQuantityA;
  typedef typename Utility::UnitTraits<DepUnitB>::template GetQuantityType<double>::type DepQuantityB;

  initialize<Utility::LinLin>( distribution );

  // Copy from unitless distribution to distribution type A
  Utility::UnitAwareAnalogElasticOneDDistribution<Utility::LinLin,IndepUnitA,DepUnitA>
    unit_aware_dist_a_copy = Utility::UnitAwareAnalogElasticOneDDistribution<Utility::LinLin,IndepUnitA,DepUnitA>::fromUnitlessDistribution( *Teuchos::rcp_dynamic_cast<Utility::AnalogElasticOneDDistribution<Utility::LinLin> >( distribution ) );

  // Copy from distribution type A to distribution type B
  Utility::UnitAwareAnalogElasticOneDDistribution<Utility::LinLin,IndepUnitB,DepUnitB>
    unit_aware_dist_b_copy( unit_aware_dist_a_copy );

  IndepQuantityA indep_quantity_a =
    Utility::QuantityTraits<IndepQuantityA>::initializeQuantity( 0.0 );
  InverseIndepQuantityA inv_indep_quantity_a =
    Utility::QuantityTraits<InverseIndepQuantityA>::initializeQuantity( 0.017233951046594968293 );
  DepQuantityA dep_quantity_a =
    Utility::QuantityTraits<DepQuantityA>::initializeQuantity( 1.0 );

  IndepQuantityB indep_quantity_b( indep_quantity_a );
  InverseIndepQuantityB inv_indep_quantity_b( inv_indep_quantity_a );
  DepQuantityB dep_quantity_b( dep_quantity_a );

  UTILITY_TEST_FLOATING_EQUALITY(
                           unit_aware_dist_a_copy.evaluate( indep_quantity_a ),
                           dep_quantity_a,
                           1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY(
                        unit_aware_dist_a_copy.evaluatePDF( indep_quantity_a ),
                        inv_indep_quantity_a,
                        1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY(
                           unit_aware_dist_b_copy.evaluate( indep_quantity_b ),
                           dep_quantity_b,
                           1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY(
                        unit_aware_dist_b_copy.evaluatePDF( indep_quantity_b ),
                        inv_indep_quantity_b,
                        1e-15 );

  Utility::setQuantity( indep_quantity_a, 0.5 );
  Utility::setQuantity( inv_indep_quantity_a, 0.0017233951046594968293 );
  Utility::setQuantity( dep_quantity_a, 0.1 );

  indep_quantity_b = IndepQuantityB( indep_quantity_a );
  inv_indep_quantity_b = InverseIndepQuantityB( inv_indep_quantity_a );
  dep_quantity_b = DepQuantityB( dep_quantity_a );

  UTILITY_TEST_FLOATING_EQUALITY(
                           unit_aware_dist_a_copy.evaluate( indep_quantity_a ),
                           dep_quantity_a,
                           1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY(
                        unit_aware_dist_a_copy.evaluatePDF( indep_quantity_a ),
                        inv_indep_quantity_a,
                        1e-6 );
  UTILITY_TEST_FLOATING_EQUALITY(
                           unit_aware_dist_b_copy.evaluate( indep_quantity_b ),
                           dep_quantity_b,
                           1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY(
                        unit_aware_dist_b_copy.evaluatePDF( indep_quantity_b ),
                        inv_indep_quantity_b,
                        1e-6 );
}

typedef si::energy si_energy;
typedef cgs::energy cgs_energy;
typedef si::amount si_amount;
typedef si::length si_length;
typedef cgs::length cgs_length;
typedef si::mass si_mass;
typedef cgs::mass cgs_mass;
typedef si::dimensionless si_dimensionless;
typedef cgs::dimensionless cgs_dimensionless;

TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticOneDDistribution,
                                      explicit_conversion,
                                      si_dimensionless,
                                      si_amount,
                                      cgs_dimensionless,
                                      si_amount );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticOneDDistribution,
                                      explicit_conversion,
                                      cgs_dimensionless,
                                      si_amount,
                                      si_dimensionless,
                                      si_amount );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticOneDDistribution,
                                      explicit_conversion,
                                      si_dimensionless,
                                      si_length,
                                      cgs_dimensionless,
                                      cgs_length );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticOneDDistribution,
                                      explicit_conversion,
                                      cgs_dimensionless,
                                      cgs_length,
                                      si_dimensionless,
                                      si_length );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticOneDDistribution,
                                      explicit_conversion,
                                      si_dimensionless,
                                      si_mass,
                                      cgs_dimensionless,
                                      cgs_mass );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticOneDDistribution,
                                      explicit_conversion,
                                      cgs_dimensionless,
                                      cgs_mass,
                                      si_dimensionless,
                                      si_mass );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticOneDDistribution,
                                      explicit_conversion,
                                      si_dimensionless,
                                      si_dimensionless,
                                      cgs_dimensionless,
                                      cgs_dimensionless );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticOneDDistribution,
                                      explicit_conversion,
                                      cgs_dimensionless,
                                      cgs_dimensionless,
                                      si_dimensionless,
                                      si_dimensionless );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticOneDDistribution,
                                      explicit_conversion,
                                      si_dimensionless,
                                      void,
                                      cgs_dimensionless,
                                      void );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticOneDDistribution,
                                      explicit_conversion,
                                      cgs_dimensionless,
                                      void,
                                      si_dimensionless,
                                      void );

//---------------------------------------------------------------------------//
// Custom setup
//---------------------------------------------------------------------------//
UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_BEGIN();

std::string test_dists_xml_file;

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_COMMAND_LINE_OPTIONS()
{
  clp().setOption( "test_dists_xml_file",
                   &test_dists_xml_file,
                   "Test distributions xml file name" );
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_DATA_INITIALIZATION()
{
  TEUCHOS_ADD_TYPE_CONVERTER( Utility::AnalogElasticOneDDistribution<Utility::LinLin> );
  TEUCHOS_ADD_TYPE_CONVERTER( Utility::AnalogElasticOneDDistribution<Utility::LogLin> );

  typedef Utility::UnitAwareAnalogElasticOneDDistribution<Utility::LinLin,si::dimensionless,si::amount> UnitAwareAnalogElasticLinLinDist;
  TEUCHOS_ADD_TYPE_CONVERTER( UnitAwareAnalogElasticLinLinDist );
  typedef Utility::UnitAwareAnalogElasticOneDDistribution<Utility::LogLin,si::dimensionless,si::amount> UnitAwareAnalogElasticLogLinDist;
  TEUCHOS_ADD_TYPE_CONVERTER( UnitAwareAnalogElasticLogLinDist );

  test_dists_list = Teuchos::getParametersFromXmlFile( test_dists_xml_file );

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstAnalogElasticOneDDistribution.cpp
//---------------------------------------------------------------------------//
