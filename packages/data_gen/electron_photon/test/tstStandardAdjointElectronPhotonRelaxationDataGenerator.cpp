//---------------------------------------------------------------------------//
//!
//! \file   tstStandardAdjointElectronPhotonRelaxationDataGenerator.cpp
//! \author Luke Kersting
//! \brief  Standard adjoint electron-photon-relaxation data generator unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Boost Includes
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/unordered_map.hpp>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

// FRENSIE Includes
#include "DataGen_StandardAdjointElectronPhotonRelaxationDataGenerator.hpp"
#include "Data_AdjointElectronPhotonRelaxationVolatileDataContainer.hpp"
#include "Data_ElectronPhotonRelaxationDataContainer.hpp"
#include "Utility_UnitTestHarnessExtensions.hpp"

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//

std::shared_ptr<Data::ElectronPhotonRelaxationDataContainer>
  h_epr_data_container, c_epr_data_container;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Check that a data generator can be constructed
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   basic_constructor )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  TEST_EQUALITY_CONST( generator.getAtomicNumber(), 1 );
  TEST_EQUALITY_CONST( generator.getMinPhotonEnergy(), 1e-3 );
  TEST_EQUALITY_CONST( generator.getMaxPhotonEnergy(), 20.0 );
  TEST_EQUALITY_CONST( generator.getMinElectronEnergy(), 1e-5 );
  TEST_EQUALITY_CONST( generator.getMaxElectronEnergy(), 1e5 );
  TEST_EQUALITY_CONST( generator.getDefaultGridConvergenceTolerance(), 1e-3 );
  TEST_EQUALITY_CONST( generator.getDefaultGridAbsoluteDifferenceTolerance(),
                       1e-12 );
  TEST_EQUALITY_CONST( generator.getDefaultGridDistanceTolerance(), 1e-14 );
  TEST_EQUALITY_CONST( generator.getAdjointPairProductionEnergyDistNormConstEvaluationTolerance(),
                       1e-3 );
  TEST_EQUALITY_CONST( generator.getAdjointPairProductionEnergyDistNormConstNudgeValue(),
                       1e-6 );
  TEST_EQUALITY_CONST( generator.getAdjointTripletProductionEnergyDistNormConstEvaluationTolerance(),
                       1e-3 );
  TEST_EQUALITY_CONST( generator.getAdjointTripletProductionEnergyDistNormConstNudgeValue(),
                       1e-6 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentMaxEnergyNudgeValue(),
                       0.2 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentEnergyToMaxEnergyNudgeValue(),
                       1e-6 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentCrossSectionEvaluationTolerance(),
                       1e-3 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentGridConvergenceTolerance(),
                       1e-3 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentGridAbsoluteDifferenceTolerance(),
                       1e-20 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentGridDistanceTolerance(),
                       1e-14 );
  TEST_EQUALITY_CONST( generator.getCutoffAngleCosine(), 1.0 );
  TEST_EQUALITY_CONST( generator.getNumberOfMomentPreservingAngles(), 0 );
  TEST_EQUALITY_CONST( generator.getAdjointBremsstrahlungCrossSectionEvaluationTolerance(),
                       1e-3 );
}

//---------------------------------------------------------------------------//
// Check that a data generator can be constructed
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   constructor )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( c_epr_data_container, 1e-3, 20.0, 1e-5, 1e5 );

  TEST_EQUALITY_CONST( generator.getAtomicNumber(), 6 );
  TEST_EQUALITY_CONST( generator.getMinPhotonEnergy(), 1e-3 );
  TEST_EQUALITY_CONST( generator.getMaxPhotonEnergy(), 20.0 );
  TEST_EQUALITY_CONST( generator.getMinElectronEnergy(), 1e-5 );
  TEST_EQUALITY_CONST( generator.getMaxElectronEnergy(), 1e5 );
  TEST_EQUALITY_CONST( generator.getDefaultGridConvergenceTolerance(), 1e-3 );
  TEST_EQUALITY_CONST( generator.getDefaultGridAbsoluteDifferenceTolerance(),
                       1e-12 );
  TEST_EQUALITY_CONST( generator.getDefaultGridDistanceTolerance(), 1e-14 );
  TEST_EQUALITY_CONST( generator.getAdjointPairProductionEnergyDistNormConstEvaluationTolerance(),
                       1e-3 );
  TEST_EQUALITY_CONST( generator.getAdjointPairProductionEnergyDistNormConstNudgeValue(),
                       1e-6 );
  TEST_EQUALITY_CONST( generator.getAdjointTripletProductionEnergyDistNormConstEvaluationTolerance(),
                       1e-3 );
  TEST_EQUALITY_CONST( generator.getAdjointTripletProductionEnergyDistNormConstNudgeValue(),
                       1e-6 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentMaxEnergyNudgeValue(),
                       0.2 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentEnergyToMaxEnergyNudgeValue(),
                       1e-6 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentCrossSectionEvaluationTolerance(),
                       1e-3 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentGridConvergenceTolerance(),
                       1e-3 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentGridAbsoluteDifferenceTolerance(),
                       1e-20 );
  TEST_EQUALITY_CONST( generator.getAdjointIncoherentGridDistanceTolerance(),
                       1e-14 );
  TEST_EQUALITY_CONST( generator.getCutoffAngleCosine(), 1.0 );
  TEST_EQUALITY_CONST( generator.getNumberOfMomentPreservingAngles(), 0 );
  TEST_EQUALITY_CONST( generator.getAdjointBremsstrahlungCrossSectionEvaluationTolerance(),
                       1e-3 );
}

//---------------------------------------------------------------------------//
// Check that the default grid convergence tolerance can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setDefaultGridConvergenceTolerance )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setDefaultGridConvergenceTolerance( 1e-5 );

  TEST_EQUALITY_CONST( generator.getDefaultGridConvergenceTolerance(), 1e-5 );
}

//---------------------------------------------------------------------------//
// Check that the default grid absolute difference tolerance can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setDefaultGridAbsoluteDifferenceTolerance )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setDefaultGridAbsoluteDifferenceTolerance( 1e-40 );

  TEST_EQUALITY_CONST( generator.getDefaultGridAbsoluteDifferenceTolerance(),
                       1e-40 );
}

//---------------------------------------------------------------------------//
// Check that the default grid distance tolerance can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setDefaultGridDistanceTolerance )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setDefaultGridDistanceTolerance( 1e-30 );

  TEST_EQUALITY_CONST( generator.getDefaultGridDistanceTolerance(), 1e-30 );
}

//---------------------------------------------------------------------------//
// Check that the adjoint pair production energy dist. norm constant evaluation
// tolerance can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setAdjointPairProductionEnergyDistNormConstEvaluationTolerance )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setAdjointPairProductionEnergyDistNormConstEvaluationTolerance( 1e-5 );

  TEST_EQUALITY_CONST( generator.getAdjointPairProductionEnergyDistNormConstEvaluationTolerance(),
                       1e-5 );
}

//---------------------------------------------------------------------------//
// Check that the adjoint pair production energy dist. norm constant nudge
// value can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setAdjointPairProductionEnergyDistNormConstNudgeValue )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setAdjointPairProductionEnergyDistNormConstNudgeValue( 1e-4 );

  TEST_EQUALITY_CONST( generator.getAdjointPairProductionEnergyDistNormConstNudgeValue(),
                       1e-4 );
}

//---------------------------------------------------------------------------//
// Check that the adjoint triplet production energy dist. norm constant
// evaluation tolerance can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setAdjointTripletProductionEnergyDistNormConstEvaluationTolerance )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setAdjointTripletProductionEnergyDistNormConstEvaluationTolerance( 1e-4 );

  TEST_EQUALITY_CONST( generator.getAdjointTripletProductionEnergyDistNormConstEvaluationTolerance(),
                       1e-4 );
}

//---------------------------------------------------------------------------//
// Check that the adjoint triplet production energy dist. nomr const
// nudge value can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setAdjointTripletProductionEnergyDistNormConstNudgeValue )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setAdjointTripletProductionEnergyDistNormConstNudgeValue( 1e-2 );

  TEST_EQUALITY_CONST( generator.getAdjointTripletProductionEnergyDistNormConstNudgeValue(),
                       1e-2 );
}

//---------------------------------------------------------------------------//
// Check that the adjoint incoherent max energy nudge value can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setAdjointIncoherentMaxEnergyNudgeValue )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setAdjointIncoherentMaxEnergyNudgeValue( 0.5 );

  TEST_EQUALITY_CONST( generator.getAdjointIncoherentMaxEnergyNudgeValue(),
                       0.5 );
}

//---------------------------------------------------------------------------//
// Check that the adjoint incoherent energy to max energy nudge value can
// be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setAdjointIncoherentEnergyToMaxEnergyNudgeValue )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setAdjointIncoherentEnergyToMaxEnergyNudgeValue( 1e-3 );

  TEST_EQUALITY_CONST( generator.getAdjointIncoherentEnergyToMaxEnergyNudgeValue(),
                       1e-3 );
}

//---------------------------------------------------------------------------//
// Check that the adjoint incoherent cross section evaluation tolerance can
// be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setAdjointIncoherentCrossSectionEvaluationTolerance )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setAdjointIncoherentCrossSectionEvaluationTolerance( 1e-4 );

  TEST_EQUALITY_CONST( generator.getAdjointIncoherentCrossSectionEvaluationTolerance(),
                       1e-4 );
}

//---------------------------------------------------------------------------//
// Check that the adjoint incoherent grid convergence tolerance can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setAdjointIncoherentGridConvergenceTolerance )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setAdjointIncoherentGridConvergenceTolerance( 1e-5 );

  TEST_EQUALITY_CONST( generator.getAdjointIncoherentGridConvergenceTolerance(),
                       1e-5 );
}

//---------------------------------------------------------------------------//
// Check that the adjoint incoherent grid absolute difference tolerance can
// be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setAdjointIncoherentGridAbsoluteDifferenceTolerance )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setAdjointIncoherentGridAbsoluteDifferenceTolerance( 1e-50 );

  TEST_EQUALITY_CONST( generator.getAdjointIncoherentGridAbsoluteDifferenceTolerance(),
                       1e-50 );
}

//---------------------------------------------------------------------------//
// Check that the adjoint incoherent grid distance tolerance can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setAdjointIncoherentGridDistanceTolerance )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setAdjointIncoherentGridDistanceTolerance( 1e-40 );

  TEST_EQUALITY_CONST( generator.getAdjointIncoherentGridDistanceTolerance(),
                       1e-40 );
}

//---------------------------------------------------------------------------//
// Check that the cutoff angle cosine can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setCutoffAngleCosine )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setCutoffAngleCosine( 0.99 );

  TEST_EQUALITY_CONST( generator.getCutoffAngleCosine(), 0.99 );
}

//---------------------------------------------------------------------------//
// Check that the number of moment preserving angles can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setNumberOfMomentPreservingAngles )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setNumberOfMomentPreservingAngles( 4 );

  TEST_EQUALITY_CONST( generator.getNumberOfMomentPreservingAngles(), 4 );
}

//---------------------------------------------------------------------------//
// Check that the adjoint bremsstrahlung evaluation tolerance can be set
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
                   setAdjointBremsstrahlungCrossSectionEvaluationTolerance )
{
  DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
    generator( h_epr_data_container );

  generator.setAdjointBremsstrahlungCrossSectionEvaluationTolerance( 1e-4 );

  TEST_EQUALITY_CONST( generator.getAdjointBremsstrahlungCrossSectionEvaluationTolerance(),
                       1e-4 );
}

//---------------------------------------------------------------------------//
// Check that a data container can be populated
TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
		   populateEPRDataContainer_h )
{
  Data::AdjointElectronPhotonRelaxationVolatileDataContainer data_container;

  {
    DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator
      generator( h_epr_data_container );
    
    generator.setDefaultGridConvergenceTolerance( 1e-3 );
    generator.setDefaultGridAbsoluteDifferenceTolerance( 1e-42 );
    generator.setDefaultGridDistanceTolerance( 1e-15 );
    generator.setAdjointPairProductionEnergyDistNormConstEvaluationTolerance( 1e-3 );
    generator.setAdjointPairProductionEnergyDistNormConstNudgeValue( 1e-6 );
    generator.setAdjointTripletProductionEnergyDistNormConstEvaluationTolerance( 1e-3 );
    generator.setAdjointTripletProductionEnergyDistNormConstNudgeValue( 1e-6 );
    generator.setAdjointIncoherentMaxEnergyNudgeValue( 0.2 );
    generator.setAdjointIncoherentEnergyToMaxEnergyNudgeValue( 1e-6 );
    generator.setAdjointIncoherentCrossSectionEvaluationTolerance( 1e-3 );
    generator.setAdjointIncoherentGridConvergenceTolerance( 0.5 );
    generator.setAdjointIncoherentGridAbsoluteDifferenceTolerance( 1e-42 );
    generator.setAdjointIncoherentGridDistanceTolerance( 1e-15 );
    generator.setCutoffAngleCosine( 0.9 );
    generator.setNumberOfMomentPreservingAngles( 1 );
    generator.setAdjointBremsstrahlungCrossSectionEvaluationTolerance( 1e-3 );

    generator.populateEPRDataContainer( data_container );
  }

  // Check the table settings data
  TEST_EQUALITY_CONST( data_container.getAtomicNumber(), 1 );
  TEST_EQUALITY_CONST( data_container.getMinPhotonEnergy(), 0.001 );
  TEST_EQUALITY_CONST( data_container.getMaxPhotonEnergy(), 20.0 );
  TEST_EQUALITY_CONST( data_container.getMinElectronEnergy(), 1.0e-5 );
  TEST_EQUALITY_CONST( data_container.getMaxElectronEnergy(), 1.0e+5 );
  TEST_EQUALITY_CONST( data_container.getGridConvergenceTolerance(), 1e-3 );
  TEST_EQUALITY_CONST( data_container.getGridAbsoluteDifferenceTolerance(),
                       1e-42 );
  TEST_EQUALITY_CONST( data_container.getGridDistanceTolerance(), 1e-15 );
  TEST_EQUALITY_CONST( data_container.getAdjointPairProductionEnergyDistNormConstantEvaluationTolerance(),
                       1e-3 );
  TEST_EQUALITY_CONST( data_container.getAdjointPairProductionEnergyDistNormConstantNudgeValue(),
                       1e-6 );
  TEST_EQUALITY_CONST( data_container.getAdjointTripletProductionEnergyDistNormConstantEvaluationTolerance(),
                       1e-3 );
  TEST_EQUALITY_CONST( data_container.getAdjointTripletProductionEnergyDistNormConstantNudgeValue(),
                       1e-6 );
  TEST_EQUALITY_CONST( data_container.getAdjointIncoherentMaxEnergyNudgeValue(),
                       0.2 );
  TEST_EQUALITY_CONST( data_container.getAdjointIncoherentEnergyToMaxEnergyNudgeValue(),
                       1e-6 );
  TEST_EQUALITY_CONST( data_container.getAdjointIncoherentEvaluationTolerance(),
                       1e-3 );
  TEST_EQUALITY_CONST( data_container.getAdjointIncoherentGridConvergenceTolerance(),
                       0.5 );
  TEST_EQUALITY_CONST( data_container.getAdjointIncoherentGridAbsoluteDifferenceTolerance(),
                       1e-42 );
  TEST_EQUALITY_CONST( data_container.getAdjointIncoherentGridDistanceTolerance(),
                       1e-15 );
  TEST_EQUALITY_CONST( data_container.getCutoffAngleCosine(), 0.9 );
  TEST_EQUALITY_CONST( data_container.getNumberOfAdjointMomentPreservingAngles(), 1 );

  // Check the relaxation data
  TEST_EQUALITY_CONST( data_container.getSubshells().size(), 1 );
  TEST_ASSERT( data_container.getSubshells().count( 1 ) );
  TEST_EQUALITY_CONST( data_container.getSubshellOccupancy( 1 ), 1 );
  TEST_EQUALITY_CONST( data_container.getSubshellBindingEnergy( 1 ),
		       1.361000000000E-05 );

  // Check the photon data

  // Check the electron data
  // std::vector<double> energy_grid = data_container.getAdjointElectronEnergyGrid();
//   TEST_EQUALITY_CONST( energy_grid.front(), 1.0e-5 );
//   TEST_EQUALITY_CONST( energy_grid.back(), 1.0e+5 );
//   TEST_EQUALITY_CONST( energy_grid.size(), 852 );

//   // Check the elastic data
//   unsigned threshold =
//     data_container.getAdjointCutoffElasticCrossSectionThresholdEnergyIndex();

//   TEST_EQUALITY_CONST( threshold, 0 );

//   std::vector<double> cross_section =
//     data_container.getAdjointCutoffElasticCrossSection();

//   TEST_EQUALITY_CONST( cross_section.front(), 2.74896e+8 );
//   TEST_FLOATING_EQUALITY( cross_section.back(), 1.31176e-5, 1e-15 );
//   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

//   threshold =
//     data_container.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex();

//   TEST_EQUALITY_CONST( threshold, 263 );

//   cross_section =
//     data_container.getAdjointScreenedRutherfordElasticCrossSection();

//   TEST_EQUALITY_CONST( cross_section.front(), 3.60034169494247180 );
//   TEST_EQUALITY_CONST( cross_section.back(), 1.29871e+4-1.31176e-5 );
//   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

//   std::vector<double> angular_grid =
//     data_container.getAdjointElasticAngularEnergyGrid();

//   TEST_EQUALITY_CONST( angular_grid.front(), 1.0e-5 );
//   TEST_EQUALITY_CONST( angular_grid.back(), 1.0e+5 );
//   TEST_EQUALITY_CONST( angular_grid.size(), 16 );

//   std::vector<double> elastic_angles =
//     data_container.getAdjointCutoffElasticAngles(1.0e-5);

//   TEST_EQUALITY_CONST( elastic_angles.front(), -1.0 );
//   TEST_EQUALITY_CONST( elastic_angles.back(), 0.999999 );
//   TEST_EQUALITY_CONST( elastic_angles.size(), 2 );

//   elastic_angles =
//     data_container.getAdjointCutoffElasticAngles(1.0e+5);

//   TEST_EQUALITY_CONST( elastic_angles.front(), -1.0 );
//   TEST_EQUALITY_CONST( elastic_angles.back(), 0.999999 );
//   TEST_EQUALITY_CONST( elastic_angles.size(), 96 );

//   std::vector<double> elastic_pdf =
//     data_container.getAdjointCutoffElasticPDF(1.0e-5);

//   TEST_EQUALITY_CONST( elastic_pdf.front(), 0.5 );
//   TEST_EQUALITY_CONST( elastic_pdf.back(), 0.5 );
//   TEST_EQUALITY_CONST( elastic_pdf.size(), 2 );

//   elastic_pdf =
//     data_container.getAdjointCutoffElasticPDF(1.0e+5);

//   TEST_EQUALITY_CONST( elastic_pdf.front(), 6.25670e-13 );
//   TEST_EQUALITY_CONST( elastic_pdf.back(), 9.86945e+5 );
//   TEST_EQUALITY_CONST( elastic_pdf.size(), 96 );

//   TEST_ASSERT( data_container.hasAdjointMomentPreservingData() );

//   std::vector<double> discrete_angles =
//     data_container.getAdjointMomentPreservingElasticDiscreteAngles( 1.0e-5 );

//   TEST_EQUALITY_CONST( discrete_angles.front(), 9.33333333326667125e-01 );
//   TEST_EQUALITY_CONST( discrete_angles.back(), 9.33333333326667125e-01 );
//   TEST_EQUALITY_CONST( discrete_angles.size(), 1 );

//   discrete_angles =
//     data_container.getAdjointMomentPreservingElasticDiscreteAngles( 1.0e+5 );

//   TEST_EQUALITY_CONST( discrete_angles.front(), 9.96847743255378838e-01 );
//   TEST_EQUALITY_CONST( discrete_angles.back(), 9.96847743255378838e-01 );
//   TEST_EQUALITY_CONST( discrete_angles.size(), 1 );

//   std::vector<double> discrete_weights =
//     data_container.getAdjointMomentPreservingElasticWeights( 1.0e-5 );

//   TEST_EQUALITY_CONST( discrete_weights.front(), 1.0 );
//   TEST_EQUALITY_CONST( discrete_weights.back(), 1.0 );
//   TEST_EQUALITY_CONST( discrete_weights.size(), 1 );

//   discrete_weights =
//     data_container.getAdjointMomentPreservingElasticWeights( 1.0e+5 );

//   TEST_EQUALITY_CONST( discrete_weights.front(), 1.0 );
//   TEST_EQUALITY_CONST( discrete_weights.back(), 1.0 );
//   TEST_EQUALITY_CONST( discrete_weights.size(), 1 );

//   threshold =
//     data_container.getAdjointMomentPreservingCrossSectionThresholdEnergyIndex();

//   TEST_EQUALITY_CONST( threshold, 0 );

//   cross_section =
//     data_container.getAdjointMomentPreservingCrossSection();

//   TEST_FLOATING_EQUALITY( cross_section.front(), 1.0308605152240909636E+07, 1e-15 );
//   TEST_FLOATING_EQUALITY( cross_section.back(), 1.2931601408114005462e-07, 1e-15 );
//   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

//   // Check the bremsstrahlung data
//   threshold =
//     data_container.getAdjointBremsstrahlungCrossSectionThresholdEnergyIndex();

//   TEST_EQUALITY_CONST( threshold, 0 );

//   cross_section =
//     data_container.getAdjointBremsstrahlungCrossSection();

//   TEST_EQUALITY_CONST( cross_section.front(), 4.34999554761655176e+01 );
//   TEST_EQUALITY_CONST( cross_section.back(), 0.0 );
//   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

//   std::vector<double> bremsstrahlung_energy_grid =
//     data_container.getAdjointBremsstrahlungEnergyGrid();

//   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.front(), 1.00000e-5 );
//   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.back(), 1.00000e+5 );
//   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.size(), 11 );

//   std::vector<double> bremsstrahlung_photon_energy =
//     data_container.getAdjointBremsstrahlungPhotonEnergy( 1.00000e-5 );

//   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.front(), 1.00000e-7 );
//   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.back(), 1.00000e-5 );
//   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.size(), 17 );

//   bremsstrahlung_photon_energy =
//     data_container.getAdjointBremsstrahlungPhotonEnergy( 1.00000e+5 );

//   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.front(), 1.00000e-7 );
//   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.back(), 1.00000e+5 );
//   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.size(), 2 );

//   std::vector<double> bremsstrahlung_photon_pdf =
//     data_container.getAdjointBremsstrahlungPhotonPDF( 1.00000e-5 );

//   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.front(), 1.46394848379322421e+06 );
//   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.back(), 1.46824029287934118e+04 );
//   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.size(), 17 );

//   bremsstrahlung_photon_pdf =
//     data_container.getAdjointBremsstrahlungPhotonPDF( 1.00000e+5 );

//   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.front(),  0.0 );
//   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.back(),  0.0 );
//   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.size(), 2 );

// /*
//   // Check the electroionization data
//   threshold =
//     data_container.getAdjointElectroionizationCrossSectionThresholdEnergyIndex( 1u );

//   TEST_EQUALITY_CONST( threshold, 7 );

//   cross_section =
//     data_container.getAdjointElectroionizationCrossSection( 1u );

// //  TEST_EQUALITY_CONST( cross_section.front(), 1.26041968911917554e+06 );
// //! \todo double check what the front cross section should be
//   TEST_EQUALITY_CONST( cross_section.front(), 1.26041968911917577e+06 );
//   TEST_EQUALITY_CONST( cross_section.back(), 8.28924e+4 );
//   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

//   std::vector<double> electroionization_energy_grid =
//     data_container.getAdjointElectroionizationEnergyGrid( 1u );

//   TEST_EQUALITY_CONST( electroionization_energy_grid.front(), 1.36100e-5 );
//   TEST_EQUALITY_CONST( electroionization_energy_grid.back(), 1.00000e+5 );
//   TEST_EQUALITY_CONST( electroionization_energy_grid.size(), 8 );

//   std::vector<double> electroionization_recoil_energy =
//     data_container.getAdjointElectroionizationRecoilEnergy( 1u, 1.36100e-5 );

//   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 2.79866e-9 );
//   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 2.79866e-8 );
//   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 2 );

//   electroionization_recoil_energy =
//     data_container.getAdjointElectroionizationRecoilEnergy( 1u, 1.00000e+5 );

//   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 1.00000e-7 );
//   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 5.00000e+4 );
//   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 147 );

//   std::vector<double> electroionization_recoil_pdf =
//     data_container.getAdjointElectroionizationRecoilPDF( 1u, 1.36100e-5 );

//   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 3.97015e+7 );
//   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 3.97015e+7 );
//   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 2 );

//   electroionization_recoil_pdf =
//     data_container.getAdjointElectroionizationRecoilPDF( 1u, 1.00000e+5 );

//   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 1.61897e+5 );
//   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 2.77550e-15 );
//   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 147 );
// */
//   // Check the atomic excitation data
//   threshold =
//     data_container.getAdjointAtomicExcitationCrossSectionThresholdEnergyIndex();

//   TEST_EQUALITY_CONST( threshold, 0 );

//   cross_section =
//     data_container.getAdjointAtomicExcitationCrossSection();

//   TEST_EQUALITY_CONST( cross_section.front(), 6.12229969785753563e+07 );
//   TEST_EQUALITY_CONST( cross_section.back(), 0.0 );
//   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

//   std::vector<double> atomic_excitation_energy_grid =
//     data_container.getAdjointAtomicExcitationEnergyGrid();

//   TEST_FLOATING_EQUALITY( atomic_excitation_energy_grid.front(), 8.7e-8, 1e-13 );
//   TEST_FLOATING_EQUALITY( atomic_excitation_energy_grid.back(), 1.0e+5 - 2.10777e-5, 1e-15 );
//   TEST_EQUALITY_CONST( atomic_excitation_energy_grid.size(), 169 );

//   std::vector<double> atomic_excitation_energy_gain =
//     data_container.getAdjointAtomicExcitationEnergyGain();

//   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.front(), 1.3697e-5 );
//   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.back(), 2.10777e-5 );
//   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.size(), 169 );

  data_container.exportData( "test_h_epr.xml",
			     Utility::ArchivableObject::XML_ARCHIVE );
}

//---------------------------------------------------------------------------//
// Check that a data container can be populated
// TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
// 		   repopulateAdjointMomentPreservingData_h )
// {
//   Data::AdjointElectronPhotonRelaxationVolatileDataContainer
//     data_container( "test_h_epr.xml",
//                     Utility::ArchivableObject::XML_ARCHIVE );

//   DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator::repopulateAdjointMomentPreservingData(
//     data_container,
//     0.9,
//     2 );

//   // Check the table settings data
//   TEST_EQUALITY_CONST( data_container.getAtomicNumber(), 1 );
//   TEST_EQUALITY_CONST( data_container.getMinPhotonEnergy(), 0.001 );
//   TEST_EQUALITY_CONST( data_container.getMaxPhotonEnergy(), 20.0 );
//   TEST_EQUALITY_CONST( data_container.getMinElectronEnergy(), 1.0e-5 );
//   TEST_EQUALITY_CONST( data_container.getMaxElectronEnergy(), 1.0e+5 );
//   TEST_EQUALITY_CONST( data_container.getCutoffAngleCosine(), 0.9 );
//   TEST_EQUALITY_CONST( data_container.getNumberOfAdjointMomentPreservingAngles(), 2 );
//   TEST_EQUALITY_CONST( data_container.getGridConvergenceTolerance(), 0.001 );
//   TEST_EQUALITY_CONST(
//     data_container.getGridAbsoluteDifferenceTolerance(), 1e-42 );
//   TEST_EQUALITY_CONST( data_container.getGridDistanceTolerance(), 1e-15 );

//   // Check the relaxation data
//   TEST_EQUALITY_CONST( data_container.getSubshells().size(), 1 );
//   TEST_ASSERT( data_container.getSubshells().count( 1 ) );
//   TEST_EQUALITY_CONST( data_container.getSubshellOccupancy( 1 ), 1 );
//   TEST_EQUALITY_CONST( data_container.getSubshellBindingEnergy( 1 ),
// 		       1.361000000000E-05 );
//   TEST_ASSERT( !data_container.hasAdjointRelaxationData() );


//   // Check the photon data

//   // Check the electron data
//   // std::vector<double> energy_grid = data_container.getAdjointElectronEnergyGrid();
// //   TEST_EQUALITY_CONST( energy_grid.front(), 1.0e-5 );
// //   TEST_EQUALITY_CONST( energy_grid.back(), 1.0e+5 );
// //   TEST_EQUALITY_CONST( energy_grid.size(), 852 );

// //   // Check the elastic data
// //   TEST_ASSERT( data_container.hasAdjointMomentPreservingData() );

// //   std::vector<double> discrete_angles =
// //     data_container.getAdjointMomentPreservingElasticDiscreteAngles( 1.0e-5 );

// //   TEST_EQUALITY_CONST( discrete_angles.front(), 9.15505102565478457e-01 );
// //   TEST_EQUALITY_CONST( discrete_angles.back(), 9.64494897399291506e-01 );
// //   TEST_EQUALITY_CONST( discrete_angles.size(), 2 );

// //   discrete_angles =
// //     data_container.getAdjointMomentPreservingElasticDiscreteAngles( 1.0e+5 );

// //   TEST_EQUALITY_CONST( discrete_angles.front(), 9.33299176569112543e-01 );
// //   TEST_EQUALITY_CONST( discrete_angles.back(), 9.99151923330158231e-01 );
// //   TEST_EQUALITY_CONST( discrete_angles.size(), 2 );

// //   std::vector<double> discrete_weights =
// //     data_container.getAdjointMomentPreservingElasticWeights( 1.0e-5 );

// //   TEST_EQUALITY_CONST( discrete_weights.front(), 4.23453445543248319e-01 );
// //   TEST_EQUALITY_CONST( discrete_weights.back(), 5.76546554456751736e-01 );
// //   TEST_EQUALITY_CONST( discrete_weights.size(), 2 );

// //   discrete_weights =
// //     data_container.getAdjointMomentPreservingElasticWeights( 1.0e+5 );

// //   TEST_EQUALITY_CONST( discrete_weights.front(), 4.60802055829827805e-04 );
// //   TEST_EQUALITY_CONST( discrete_weights.back(), 9.99539197944170232e-01 );
// //   TEST_EQUALITY_CONST( discrete_weights.size(), 2 );

// //   unsigned threshold =
// //     data_container.getAdjointMomentPreservingCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   std::vector<double> cross_section =
// //     data_container.getAdjointMomentPreservingCrossSection();

// //   TEST_FLOATING_EQUALITY( cross_section.front(), 1.2217606103336416185e+07, 1e-15 );
// //   TEST_FLOATING_EQUALITY( cross_section.back(), 4.6405644799051960388e-07, 1e-15 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

// //   threshold =
// //     data_container.getAdjointCutoffElasticCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   cross_section =
// //     data_container.getAdjointCutoffElasticCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 2.74896e+8 );
// //   TEST_FLOATING_EQUALITY( cross_section.back(), 1.31176e-5, 1e-15 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

// //   threshold =
// //     data_container.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 263 );

// //   cross_section =
// //     data_container.getAdjointScreenedRutherfordElasticCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 3.60034169494247180 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 1.29871e+4-1.31176e-5 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

// //   std::vector<double> angular_grid =
// //     data_container.getAdjointElasticAngularEnergyGrid();

// //   TEST_EQUALITY_CONST( angular_grid.front(), 1.0e-5 );
// //   TEST_EQUALITY_CONST( angular_grid.back(), 1.0e+5 );
// //   TEST_EQUALITY_CONST( angular_grid.size(), 16 );

// //   std::vector<double> elastic_angles =
// //     data_container.getAdjointCutoffElasticAngles(1.0e-5);

// //   TEST_EQUALITY_CONST( elastic_angles.front(), -1.0 );
// //   TEST_EQUALITY_CONST( elastic_angles.back(), 0.999999 );
// //   TEST_EQUALITY_CONST( elastic_angles.size(), 2 );

// //   elastic_angles =
// //     data_container.getAdjointCutoffElasticAngles(1.0e+5);

// //   TEST_EQUALITY_CONST( elastic_angles.front(), -1.0 );
// //   TEST_EQUALITY_CONST( elastic_angles.back(), 0.999999 );
// //   TEST_EQUALITY_CONST( elastic_angles.size(), 96 );

// //   std::vector<double> elastic_pdf =
// //     data_container.getAdjointCutoffElasticPDF(1.0e-5);

// //   TEST_EQUALITY_CONST( elastic_pdf.front(), 0.5 );
// //   TEST_EQUALITY_CONST( elastic_pdf.back(), 0.5 );
// //   TEST_EQUALITY_CONST( elastic_pdf.size(), 2 );

// //   elastic_pdf =
// //     data_container.getAdjointCutoffElasticPDF(1.0e+5);

// //   TEST_EQUALITY_CONST( elastic_pdf.front(), 6.25670e-13 );
// //   TEST_EQUALITY_CONST( elastic_pdf.back(), 9.86945e+5 );
// //   TEST_EQUALITY_CONST( elastic_pdf.size(), 96 );
// // /*
// //   // Check the electroionization data
// //   threshold =
// //     data_container.getAdjointElectroionizationCrossSectionThresholdEnergyIndex( 1u );

// //   TEST_EQUALITY_CONST( threshold, 7 );

// //   cross_section =
// //     data_container.getAdjointElectroionizationCrossSection( 1u );

// // //  TEST_EQUALITY_CONST( cross_section.front(), 1.26041968911917554e+06 );
// // //! \todo double check what the front cross section should be
// //   TEST_EQUALITY_CONST( cross_section.front(), 1.26041968911917577e+06 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 8.28924e+4 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

// //   std::vector<double> electroionization_energy_grid =
// //     data_container.getAdjointElectroionizationEnergyGrid( 1u );

// //   TEST_EQUALITY_CONST( electroionization_energy_grid.front(), 1.36100e-5 );
// //   TEST_EQUALITY_CONST( electroionization_energy_grid.back(), 1.00000e+5 );
// //   TEST_EQUALITY_CONST( electroionization_energy_grid.size(), 8 );

// //   std::vector<double> electroionization_recoil_energy =
// //     data_container.getAdjointElectroionizationRecoilEnergy( 1u, 1.36100e-5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 2.79866e-9 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 2.79866e-8 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 2 );

// //   electroionization_recoil_energy =
// //     data_container.getAdjointElectroionizationRecoilEnergy( 1u, 1.00000e+5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 5.00000e+4 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 147 );

// //   std::vector<double> electroionization_recoil_pdf =
// //     data_container.getAdjointElectroionizationRecoilPDF( 1u, 1.36100e-5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 3.97015e+7 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 3.97015e+7 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 2 );

// //   electroionization_recoil_pdf =
// //     data_container.getAdjointElectroionizationRecoilPDF( 1u, 1.00000e+5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 1.61897e+5 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 2.77550e-15 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 147 );
// // */
// //   // Check the bremsstrahlung data
// //   threshold =
// //     data_container.getAdjointBremsstrahlungCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   cross_section =
// //     data_container.getAdjointBremsstrahlungCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 4.34999554761655176e+01 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 0.0 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

// //   std::vector<double> bremsstrahlung_energy_grid =
// //     data_container.getAdjointBremsstrahlungEnergyGrid();

// //   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.front(), 1.00000e-5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.back(), 1.00000e+5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.size(), 11 );

// //   std::vector<double> bremsstrahlung_photon_energy =
// //     data_container.getAdjointBremsstrahlungPhotonEnergy( 1.00000e-5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.front(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.back(), 1.00000e-5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.size(), 17 );

// //   bremsstrahlung_photon_energy =
// //     data_container.getAdjointBremsstrahlungPhotonEnergy( 1.00000e+5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.front(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.back(), 1.00000e+5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.size(), 2 );

// //   std::vector<double> bremsstrahlung_photon_pdf =
// //     data_container.getAdjointBremsstrahlungPhotonPDF( 1.0e-5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.front(), 1.46394848379322421e+06 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.back(), 1.46824029287934118e+04 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.size(), 17 );

// //   bremsstrahlung_photon_pdf =
// //     data_container.getAdjointBremsstrahlungPhotonPDF( 1.0e+5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.front(),  0.0 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.back(),  0.0 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.size(), 2 );

// //   // Check the atomic excitation data
// //   threshold =
// //     data_container.getAdjointAtomicExcitationCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   cross_section =
// //     data_container.getAdjointAtomicExcitationCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 6.12229969785753563e+07 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 0.0 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 852-threshold );

// //   std::vector<double> atomic_excitation_energy_grid =
// //     data_container.getAdjointAtomicExcitationEnergyGrid();

// //   TEST_FLOATING_EQUALITY( atomic_excitation_energy_grid.front(), 8.7e-8, 1e-13 );
// //   TEST_FLOATING_EQUALITY( atomic_excitation_energy_grid.back(), 1.0e+5 - 2.10777e-5, 1e-15 );
// //   TEST_EQUALITY_CONST( atomic_excitation_energy_grid.size(), 169 );

// //   std::vector<double> atomic_excitation_energy_gain =
// //     data_container.getAdjointAtomicExcitationEnergyGain();

// //   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.front(), 1.3697e-5 );
// //   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.back(), 2.10777e-5 );
// //   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.size(), 169 );

// //   data_container.exportData( "test_h_epr.xml",
// // 			     Utility::ArchivableObject::XML_ARCHIVE );
// }

//---------------------------------------------------------------------------//
// Check that a data container can be populated
// TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
// 		   populateEPRDataContainer_c )
// {
//   data_generator_c.reset(
//        new DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator(
//                    c_epr_data_container, 1.0, 0, 0.001, 0.001, 1e-32, 1e-16) );

//   Data::AdjointElectronPhotonRelaxationVolatileDataContainer data_container;

//   data_generator_c->populateEPRDataContainer( data_container );

//   // Check the table settings data
//   TEST_EQUALITY_CONST( data_container.getAtomicNumber(), 6 );
//   TEST_EQUALITY_CONST( data_container.getMinPhotonEnergy(), 0.001 );
//   TEST_EQUALITY_CONST( data_container.getMaxPhotonEnergy(), 20.0 );
//   TEST_EQUALITY_CONST( data_container.getMinElectronEnergy(), 1.0e-5 );
//   TEST_EQUALITY_CONST( data_container.getMaxElectronEnergy(), 1.0e+5 );
//   TEST_EQUALITY_CONST( data_container.getCutoffAngleCosine(), 1.0 );
//   TEST_EQUALITY_CONST( data_container.getNumberOfAdjointMomentPreservingAngles(), 0 );
//   TEST_EQUALITY_CONST( data_container.getGridConvergenceTolerance(), 0.001 );
//   TEST_EQUALITY_CONST(
//     data_container.getGridAbsoluteDifferenceTolerance(), 1e-32 );
//   TEST_EQUALITY_CONST( data_container.getGridDistanceTolerance(), 1e-16 );

//   // Check the relaxation data
//   TEST_EQUALITY_CONST( data_container.getSubshells().size(), 4 );
//   TEST_ASSERT( data_container.getSubshells().count( 1 ) );
//   TEST_ASSERT( data_container.getSubshells().count( 2 ) );
//   TEST_ASSERT( data_container.getSubshells().count( 3 ) );
//   TEST_ASSERT( data_container.getSubshells().count( 4 ) );
//   TEST_EQUALITY_CONST( data_container.getSubshellOccupancy( 1 ), 2 );
//   TEST_EQUALITY_CONST( data_container.getSubshellOccupancy( 2 ), 2 );
//   TEST_EQUALITY_CONST( data_container.getSubshellOccupancy( 3 ), 0.67 );
//   TEST_EQUALITY_CONST( data_container.getSubshellOccupancy( 4 ), 1.33 );
//   TEST_EQUALITY_CONST( data_container.getSubshellBindingEnergy( 1 ),
// 		       2.9101e-4 );
//   TEST_EQUALITY_CONST( data_container.getSubshellBindingEnergy( 2 ),
// 		       1.7560e-5 );
//   TEST_EQUALITY_CONST( data_container.getSubshellBindingEnergy( 3 ),
// 		       8.9900e-6 );
//   TEST_EQUALITY_CONST( data_container.getSubshellBindingEnergy( 4 ),
// 		       8.9800e-6 );
//   TEST_ASSERT( !data_container.hasAdjointRelaxationData() );

//   // Check the photon data

//   // Check the electron data
//   // std::vector<double> energy_grid = data_container.getAdjointElectronEnergyGrid();
// //   TEST_EQUALITY_CONST( energy_grid.front(), 1.0e-5 );
// //   TEST_EQUALITY_CONST( energy_grid.back(), 1.0e+5 );
// //   TEST_EQUALITY_CONST( energy_grid.size(), 968 );

// //   // Check the elastic data
// //   TEST_ASSERT( !data_container.hasAdjointMomentPreservingData() );

// //   unsigned threshold =
// //     data_container.getAdjointCutoffElasticCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   std::vector<double> cross_section =
// //     data_container.getAdjointCutoffElasticCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 3.06351e+9 );
// //   TEST_FLOATING_EQUALITY( cross_section.back(), 4.72309e-4, 1e-15 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   threshold =
// //     data_container.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 292 );

// //   cross_section =
// //     data_container.getAdjointScreenedRutherfordElasticCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 1.89575980269990396e+02 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 1.40721999527690990e+05 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   std::vector<double> angular_grid =
// //     data_container.getAdjointElasticAngularEnergyGrid();

// //   TEST_EQUALITY_CONST( angular_grid.front(), 1.0e-5 );
// //   TEST_EQUALITY_CONST( angular_grid.back(), 1.0e+5 );
// //   TEST_EQUALITY_CONST( angular_grid.size(), 16 );

// //   std::vector<double> elastic_angles =
// //     data_container.getAdjointCutoffElasticAngles(1.0e-5);

// //   TEST_EQUALITY_CONST( elastic_angles.front(), -1.0 );
// //   TEST_EQUALITY_CONST( elastic_angles.back(), 0.999999 );
// //   TEST_EQUALITY_CONST( elastic_angles.size(), 2 );

// //   elastic_angles =
// //     data_container.getAdjointCutoffElasticAngles(1.0e+5);

// //   TEST_EQUALITY_CONST( elastic_angles.front(), -1.0 );
// //   TEST_EQUALITY_CONST( elastic_angles.back(), 0.999999 );
// //   TEST_EQUALITY_CONST( elastic_angles.size(), 96 );

// //   std::vector<double> elastic_pdf =
// //     data_container.getAdjointCutoffElasticPDF(1.0e-5);

// //   TEST_EQUALITY_CONST( elastic_pdf.front(), 0.5 );
// //   TEST_EQUALITY_CONST( elastic_pdf.back(), 0.5 );
// //   TEST_EQUALITY_CONST( elastic_pdf.size(), 2 );

// //   elastic_pdf =
// //     data_container.getAdjointCutoffElasticPDF(1.0e+5);

// //   TEST_EQUALITY_CONST( elastic_pdf.front(), 1.693970E-11 );
// //   TEST_EQUALITY_CONST( elastic_pdf.back(), 9.868670E+05 );
// //   TEST_EQUALITY_CONST( elastic_pdf.size(), 96 );

// //   // Check the bremsstrahlung data
// //   threshold =
// //     data_container.getAdjointBremsstrahlungCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   cross_section =
// //     data_container.getAdjointBremsstrahlungCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 9.13216405736755519e+02 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 0.0 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   std::vector<double> bremsstrahlung_energy_grid =
// //     data_container.getAdjointBremsstrahlungEnergyGrid();

// //   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.front(), 1.00000e-5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.back(), 1.00000e+5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.size(), 10 );

// //   std::vector<double> bremsstrahlung_photon_energy =
// //     data_container.getAdjointBremsstrahlungPhotonEnergy( 1.00000e-5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.front(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.back(), 1.00000e-5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.size(), 17 );

// //   bremsstrahlung_photon_energy =
// //     data_container.getAdjointBremsstrahlungPhotonEnergy( 1.00000e+5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.front(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.back(), 1.00000e+5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.size(), 2 );

// //   std::vector<double> bremsstrahlung_photon_pdf =
// //     data_container.getAdjointBremsstrahlungPhotonPDF( 1.00000e-5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.front(),1.40925689457137836e+06 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.back(), 1.43093090368334142e+04 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.size(), 17 );

// //   bremsstrahlung_photon_pdf =
// //     data_container.getAdjointBremsstrahlungPhotonPDF( 1.00000e+5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.front(), 0.0 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.back(),  0.0 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.size(), 2 );
// // /*
// //   // Check the electroionization data
// //   threshold =
// //     data_container.getAdjointElectroionizationCrossSectionThresholdEnergyIndex( 1u );

// //   TEST_EQUALITY_CONST( threshold, 85 );

// //   cross_section =
// //     data_container.getAdjointElectroionizationCrossSection( 1u );

// //   TEST_EQUALITY_CONST( cross_section.front(), 1.10703713515228592e+04 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 1.338050E+04 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   std::vector<double> electroionization_energy_grid =
// //     data_container.getAdjointElectroionizationEnergyGrid( 1u );

// //   TEST_EQUALITY_CONST( electroionization_energy_grid.front(), 2.910100E-04 );
// //   TEST_EQUALITY_CONST( electroionization_energy_grid.back(), 1.00000e+5 );
// //   TEST_EQUALITY_CONST( electroionization_energy_grid.size(), 7 );

// //   std::vector<double> electroionization_recoil_energy =
// //     data_container.getAdjointElectroionizationRecoilEnergy( 1u, 2.910100E-04 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 1.00000e-8 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 2 );

// //   electroionization_recoil_energy =
// //     data_container.getAdjointElectroionizationRecoilEnergy( 1u, 1.00000e+5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 5.00000e+4 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 128 );

// //   std::vector<double> electroionization_recoil_pdf =
// //     data_container.getAdjointElectroionizationRecoilPDF( 1u, 2.910100E-04 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 1.111110E+07 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 1.111110E+07 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 2 );

// //   electroionization_recoil_pdf =
// //     data_container.getAdjointElectroionizationRecoilPDF( 1u, 1.00000e+5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 7.358100E+03 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 3.45597E-14 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 128 );


// //   threshold =
// //     data_container.getAdjointElectroionizationCrossSectionThresholdEnergyIndex( 4u );

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   cross_section =
// //     data_container.getAdjointElectroionizationCrossSection( 4u );

// //   TEST_EQUALITY_CONST( cross_section.front(), 2.102930E+07 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 2.017010E+05 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   electroionization_energy_grid =
// //     data_container.getAdjointElectroionizationEnergyGrid( 4u );

// //   TEST_EQUALITY_CONST( electroionization_energy_grid.front(), 8.980000E-06 );
// //   TEST_EQUALITY_CONST( electroionization_energy_grid.back(), 1.00000e+5 );
// //   TEST_EQUALITY_CONST( electroionization_energy_grid.size(), 8 );

// //   electroionization_recoil_energy =
// //     data_container.getAdjointElectroionizationRecoilEnergy( 4u, 8.980000E-06 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 2.550000E-09 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 2.550000E-08 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 2 );

// //   electroionization_recoil_energy =
// //     data_container.getAdjointElectroionizationRecoilEnergy( 4u, 1.00000e+5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 5.00000e+4 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 143 );

// //   electroionization_recoil_pdf =
// //     data_container.getAdjointElectroionizationRecoilPDF( 4u, 8.980000E-06 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 4.357300E+07 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 4.357300E+07 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 2 );

// //   electroionization_recoil_pdf =
// //     data_container.getAdjointElectroionizationRecoilPDF( 4u, 1.00000e+5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 1.120930E+05 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 1.515230E-15 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 143 );
// // */
// //   // Check the atomic excitation data
// //   threshold =
// //     data_container.getAdjointAtomicExcitationCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   cross_section =
// //     data_container.getAdjointAtomicExcitationCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 5.90068663943557292e+07 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 0.0 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   std::vector<double> atomic_excitation_energy_grid =
// //     data_container.getAdjointAtomicExcitationEnergyGrid();

// //   TEST_FLOATING_EQUALITY( atomic_excitation_energy_grid.front(), 7.6731e-07, 1e-14 );
// //   TEST_FLOATING_EQUALITY( atomic_excitation_energy_grid.back(), 1.00000e+5 - 1.98154e-5, 1e-15 );
// //   TEST_EQUALITY_CONST( atomic_excitation_energy_grid.size(), 181 );

// //   std::vector<double> atomic_excitation_energy_gain =
// //     data_container.getAdjointAtomicExcitationEnergyGain();

// //   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.front(), 9.232690E-06 );
// //   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.back(), 1.981540E-05 );
// //   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.size(), 181 );

//   data_container.exportData( "test_c_epr.xml",
// 			     Utility::ArchivableObject::XML_ARCHIVE );
// }

//---------------------------------------------------------------------------//
// Check that a data container can be repopulated with moment preserving data
// TEUCHOS_UNIT_TEST( StandardAdjointElectronPhotonRelaxationDataGenerator,
//                    repopulateAdjointMomentPreservingData_c )
// {
//   Data::AdjointElectronPhotonRelaxationVolatileDataContainer
//     data_container( "test_c_epr.xml",
//                     Utility::ArchivableObject::XML_ARCHIVE );

//   DataGen::StandardAdjointElectronPhotonRelaxationDataGenerator::repopulateAdjointMomentPreservingData(
//     data_container,
//     0.9,
//     2 );

//   // Check the table settings data
//   TEST_EQUALITY_CONST( data_container.getAtomicNumber(), 6 );
//   TEST_EQUALITY_CONST( data_container.getMinPhotonEnergy(), 0.001 );
//   TEST_EQUALITY_CONST( data_container.getMaxPhotonEnergy(), 20.0 );
//   TEST_EQUALITY_CONST( data_container.getMinElectronEnergy(), 1.0e-5 );
//   TEST_EQUALITY_CONST( data_container.getMaxElectronEnergy(), 1.0e+5 );
//   TEST_EQUALITY_CONST( data_container.getCutoffAngleCosine(), 0.9 );
//   TEST_EQUALITY_CONST( data_container.getNumberOfAdjointMomentPreservingAngles(), 2 );
//   TEST_EQUALITY_CONST( data_container.getGridConvergenceTolerance(), 0.001 );
//   TEST_EQUALITY_CONST(
//     data_container.getGridAbsoluteDifferenceTolerance(), 1e-32 );
//   TEST_EQUALITY_CONST( data_container.getGridDistanceTolerance(), 1e-16 );

//   // Check the relaxation data
//   TEST_EQUALITY_CONST( data_container.getSubshells().size(), 4 );
//   TEST_ASSERT( data_container.getSubshells().count( 1 ) );
//   TEST_ASSERT( data_container.getSubshells().count( 2 ) );
//   TEST_ASSERT( data_container.getSubshells().count( 3 ) );
//   TEST_ASSERT( data_container.getSubshells().count( 4 ) );
//   TEST_EQUALITY_CONST( data_container.getSubshellOccupancy( 1 ), 2 );
//   TEST_EQUALITY_CONST( data_container.getSubshellOccupancy( 2 ), 2 );
//   TEST_EQUALITY_CONST( data_container.getSubshellOccupancy( 3 ), 0.67 );
//   TEST_EQUALITY_CONST( data_container.getSubshellOccupancy( 4 ), 1.33 );
//   TEST_EQUALITY_CONST( data_container.getSubshellBindingEnergy( 1 ),
// 		       2.9101e-4 );
//   TEST_EQUALITY_CONST( data_container.getSubshellBindingEnergy( 2 ),
// 		       1.7560e-5 );
//   TEST_EQUALITY_CONST( data_container.getSubshellBindingEnergy( 3 ),
// 		       8.9900e-6 );
//   TEST_EQUALITY_CONST( data_container.getSubshellBindingEnergy( 4 ),
// 		       8.9800e-6 );
//   TEST_ASSERT( !data_container.hasAdjointRelaxationData() );

//   // Check the photon data
  
//   // Check the electron energy grid data
//   // std::vector<double> energy_grid = data_container.getAdjointElectronEnergyGrid();
// //   TEST_EQUALITY_CONST( energy_grid.front(), 1.0e-5 );
// //   TEST_EQUALITY_CONST( energy_grid.back(), 1.0e+5 );
// //   TEST_EQUALITY_CONST( energy_grid.size(), 968 );

// //   // Check the elastic data
// //   TEST_ASSERT( data_container.hasAdjointMomentPreservingData() );

// //   std::vector<double> discrete_angles =
// //     data_container.getAdjointMomentPreservingElasticDiscreteAngles( 1.0e-5 );

// //   TEST_EQUALITY_CONST( discrete_angles.front(), 9.15505102565478457e-01 );
// //   TEST_EQUALITY_CONST( discrete_angles.back(), 9.64494897399291506e-01 );
// //   TEST_EQUALITY_CONST( discrete_angles.size(), 2 );

// //   discrete_angles =
// //     data_container.getAdjointMomentPreservingElasticDiscreteAngles( 1.0e+5 );

// //   TEST_EQUALITY_CONST( discrete_angles.front(), 9.33209295213873080e-01 );
// //   TEST_EQUALITY_CONST( discrete_angles.back(), 9.99107421216845260e-01 );
// //   TEST_EQUALITY_CONST( discrete_angles.size(), 2 );

// //   std::vector<double> discrete_weights =
// //     data_container.getAdjointMomentPreservingElasticWeights( 1.0e-5 );

// //   TEST_EQUALITY_CONST( discrete_weights.front(), 4.23453445543248319e-01 );
// //   TEST_EQUALITY_CONST( discrete_weights.back(), 5.76546554456751736e-01 );
// //   TEST_EQUALITY_CONST( discrete_weights.size(), 2 );

// //   discrete_weights =
// //     data_container.getAdjointMomentPreservingElasticWeights( 1.0e+5 );

// //   TEST_EQUALITY_CONST( discrete_weights.front(), 5.11428725729088773e-04 );
// //   TEST_EQUALITY_CONST( discrete_weights.back(), 9.99488571274270932e-01 );
// //   TEST_EQUALITY_CONST( discrete_weights.size(), 2 );

// //   unsigned threshold =
// //     data_container.getAdjointMomentPreservingCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   std::vector<double> cross_section =
// //     data_container.getAdjointMomentPreservingCrossSection();

// //   TEST_FLOATING_EQUALITY( cross_section.front(), 1.3615606801711243391e+08, 1e-15 );
// //   TEST_FLOATING_EQUALITY( cross_section.back(), 1.5258885009562140901e-05, 1e-15 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   threshold =
// //     data_container.getAdjointCutoffElasticCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   cross_section =
// //     data_container.getAdjointCutoffElasticCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 3.06351e+9 );
// //   TEST_FLOATING_EQUALITY( cross_section.back(), 4.72309e-4, 1e-15 );
// // std::cout << std::setprecision(20) << "cross_section.back() = " << cross_section.back() << std::endl;
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   threshold =
// //     data_container.getAdjointScreenedRutherfordElasticCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 292 );

// //   cross_section =
// //     data_container.getAdjointScreenedRutherfordElasticCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 1.89575980269990396e+02 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 1.40721999527690990e+05 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   std::vector<double> angular_grid =
// //     data_container.getAdjointElasticAngularEnergyGrid();

// //   TEST_EQUALITY_CONST( angular_grid.front(), 1.0e-5 );
// //   TEST_EQUALITY_CONST( angular_grid.back(), 1.0e+5 );
// //   TEST_EQUALITY_CONST( angular_grid.size(), 16 );

// //   std::vector<double> elastic_angles =
// //     data_container.getAdjointCutoffElasticAngles(1.0e-5);

// //   TEST_EQUALITY_CONST( elastic_angles.front(), -1.0 );
// //   TEST_EQUALITY_CONST( elastic_angles.back(), 0.999999 );
// //   TEST_EQUALITY_CONST( elastic_angles.size(), 2 );

// //   elastic_angles =
// //     data_container.getAdjointCutoffElasticAngles(1.0e+5);

// //   TEST_EQUALITY_CONST( elastic_angles.front(), -1.0 );
// //   TEST_EQUALITY_CONST( elastic_angles.back(), 0.999999 );
// //   TEST_EQUALITY_CONST( elastic_angles.size(), 96 );

// //   std::vector<double> elastic_pdf =
// //     data_container.getAdjointCutoffElasticPDF(1.0e-5);

// //   TEST_EQUALITY_CONST( elastic_pdf.front(), 0.5 );
// //   TEST_EQUALITY_CONST( elastic_pdf.back(), 0.5 );
// //   TEST_EQUALITY_CONST( elastic_pdf.size(), 2 );

// //   elastic_pdf =
// //     data_container.getAdjointCutoffElasticPDF(1.0e+5);

// //   TEST_EQUALITY_CONST( elastic_pdf.front(), 1.693970E-11 );
// //   TEST_EQUALITY_CONST( elastic_pdf.back(), 9.868670E+05 );
// //   TEST_EQUALITY_CONST( elastic_pdf.size(), 96 );

// //   // Check the bremsstrahlung data
// //   threshold =
// //     data_container.getAdjointBremsstrahlungCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   cross_section =
// //     data_container.getAdjointBremsstrahlungCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 9.13216405736755519e+02 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 0.0 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   std::vector<double> bremsstrahlung_energy_grid =
// //     data_container.getAdjointBremsstrahlungEnergyGrid();

// //   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.front(), 1.00000e-5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.back(), 1.00000e+5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_energy_grid.size(), 10 );

// //   std::vector<double> bremsstrahlung_photon_energy =
// //     data_container.getAdjointBremsstrahlungPhotonEnergy( 1.00000e-5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.front(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.back(), 1.00000e-5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.size(), 17 );

// //   bremsstrahlung_photon_energy =
// //     data_container.getAdjointBremsstrahlungPhotonEnergy( 1.00000e+5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.front(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.back(), 1.00000e+5 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_energy.size(), 2 );

// //   std::vector<double> bremsstrahlung_photon_pdf =
// //     data_container.getAdjointBremsstrahlungPhotonPDF( 1.00000e-5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.front(),1.40925689457137836e+06 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.back(), 1.43093090368334142e+04 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.size(), 17 );

// //   bremsstrahlung_photon_pdf =
// //     data_container.getAdjointBremsstrahlungPhotonPDF( 1.00000e+5 );

// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.front(), 0.0 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.back(),  0.0 );
// //   TEST_EQUALITY_CONST( bremsstrahlung_photon_pdf.size(), 2 );
// // /*
// //   // Check the electroionization data
// //   threshold =
// //     data_container.getAdjointElectroionizationCrossSectionThresholdEnergyIndex( 1u );

// //   TEST_EQUALITY_CONST( threshold, 85 );

// //   cross_section =
// //     data_container.getAdjointElectroionizationCrossSection( 1u );

// //   TEST_EQUALITY_CONST( cross_section.front(), 1.10703713515228592e+04 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 1.338050E+04 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   std::vector<double> electroionization_energy_grid =
// //     data_container.getAdjointElectroionizationEnergyGrid( 1u );

// //   TEST_EQUALITY_CONST( electroionization_energy_grid.front(), 2.910100E-04 );
// //   TEST_EQUALITY_CONST( electroionization_energy_grid.back(), 1.00000e+5 );
// //   TEST_EQUALITY_CONST( electroionization_energy_grid.size(), 7 );

// //   std::vector<double> electroionization_recoil_energy =
// //     data_container.getAdjointElectroionizationRecoilEnergy( 1u, 2.910100E-04 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 1.00000e-8 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 2 );

// //   electroionization_recoil_energy =
// //     data_container.getAdjointElectroionizationRecoilEnergy( 1u, 1.00000e+5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 5.00000e+4 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 128 );

// //   std::vector<double> electroionization_recoil_pdf =
// //     data_container.getAdjointElectroionizationRecoilPDF( 1u, 2.910100E-04 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 1.111110E+07 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 1.111110E+07 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 2 );

// //   electroionization_recoil_pdf =
// //     data_container.getAdjointElectroionizationRecoilPDF( 1u, 1.00000e+5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 7.358100E+03 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 3.45597E-14 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 128 );


// //   threshold =
// //     data_container.getAdjointElectroionizationCrossSectionThresholdEnergyIndex( 4u );

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   cross_section =
// //     data_container.getAdjointElectroionizationCrossSection( 4u );

// //   TEST_EQUALITY_CONST( cross_section.front(), 2.102930E+07 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 2.017010E+05 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   electroionization_energy_grid =
// //     data_container.getAdjointElectroionizationEnergyGrid( 4u );

// //   TEST_EQUALITY_CONST( electroionization_energy_grid.front(), 8.980000E-06 );
// //   TEST_EQUALITY_CONST( electroionization_energy_grid.back(), 1.00000e+5 );
// //   TEST_EQUALITY_CONST( electroionization_energy_grid.size(), 8 );

// //   electroionization_recoil_energy =
// //     data_container.getAdjointElectroionizationRecoilEnergy( 4u, 8.980000E-06 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 2.550000E-09 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 2.550000E-08 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 2 );

// //   electroionization_recoil_energy =
// //     data_container.getAdjointElectroionizationRecoilEnergy( 4u, 1.00000e+5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.front(), 1.00000e-7 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.back(), 5.00000e+4 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_energy.size(), 143 );

// //   electroionization_recoil_pdf =
// //     data_container.getAdjointElectroionizationRecoilPDF( 4u, 8.980000E-06 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 4.357300E+07 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 4.357300E+07 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 2 );

// //   electroionization_recoil_pdf =
// //     data_container.getAdjointElectroionizationRecoilPDF( 4u, 1.00000e+5 );

// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.front(), 1.120930E+05 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.back(), 1.515230E-15 );
// //   TEST_EQUALITY_CONST( electroionization_recoil_pdf.size(), 143 );
// // */
// //   // Check the atomic excitation data
// //   threshold =
// //     data_container.getAdjointAtomicExcitationCrossSectionThresholdEnergyIndex();

// //   TEST_EQUALITY_CONST( threshold, 0 );

// //   cross_section =
// //     data_container.getAdjointAtomicExcitationCrossSection();

// //   TEST_EQUALITY_CONST( cross_section.front(), 5.90068663943557292e+07 );
// //   TEST_EQUALITY_CONST( cross_section.back(), 0.0 );
// //   TEST_EQUALITY_CONST( cross_section.size(), 968-threshold );

// //   std::vector<double> atomic_excitation_energy_grid =
// //     data_container.getAdjointAtomicExcitationEnergyGrid();

// //   TEST_FLOATING_EQUALITY( atomic_excitation_energy_grid.front(), 7.6731e-07, 1e-14 );
// //   TEST_FLOATING_EQUALITY( atomic_excitation_energy_grid.back(), 1.00000e+5 - 1.98154e-5, 1e-15 );
// //   TEST_EQUALITY_CONST( atomic_excitation_energy_grid.size(), 181 );

// //   std::vector<double> atomic_excitation_energy_gain =
// //     data_container.getAdjointAtomicExcitationEnergyGain();

// //   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.front(), 9.232690E-06 );
// //   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.back(), 1.981540E-05 );
// //   TEST_EQUALITY_CONST( atomic_excitation_energy_gain.size(), 181 );

// //   data_container.exportData( "test_c_epr.xml",
// // 			     Utility::ArchivableObject::XML_ARCHIVE );
// }

//---------------------------------------------------------------------------//
// Custom main function
//---------------------------------------------------------------------------//
int main( int argc, char** argv )
{
  std::string test_h_native_file, test_c_native_file;

  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();

  clp.setOption( "test_h_native_file",
                 &test_h_native_file,
                 "Test h native file name" );
  clp.setOption( "test_c_native_file",
                 &test_c_native_file,
                 "Test c native file name" );
  
  const Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return =
    clp.parse(argc,argv);

  if ( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
    *out << "\nEnd Result: TEST FAILED" << std::endl;
    return parse_return;
  }

  // Create the native data file container for h
  h_epr_data_container.reset( new Data::ElectronPhotonRelaxationDataContainer(
                                                        test_h_native_file ) );
  c_epr_data_container.reset( new Data::ElectronPhotonRelaxationDataContainer(
                                                        test_c_native_file ) );

  // Run the unit tests
  Teuchos::GlobalMPISession mpiSession( &argc, &argv );

  const bool success = Teuchos::UnitTestRepository::runUnitTests( *out );

  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;

  clp.printFinalTimerSummary(out.ptr());

  return (success ? 0 : 1);
}

//---------------------------------------------------------------------------//
// end tstStandardAdjointElectronPhotonRelaxationDataGenerator.cpp
//---------------------------------------------------------------------------//
