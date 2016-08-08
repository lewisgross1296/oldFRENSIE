//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_ElectroionizationSubshellElectronScatteringDistributionENDLFactory.cpp
//! \author Luke Kersting
//! \brief  The electroionization subshell scattering distribution ENDL factory definition
//!
//---------------------------------------------------------------------------//

// FRENSIE Includes
#include "MonteCarlo_ElectroionizationSubshellElectronScatteringDistributionENDLFactory.hpp"
#include "Utility_TabularDistribution.hpp"
#include "Utility_ContractException.hpp"

namespace MonteCarlo{

// Create a electroionization subshell distribution
void ElectroionizationSubshellElectronScatteringDistributionENDLFactory::createElectroionizationSubshellDistribution(
	const Data::ENDLDataContainer& raw_electroionization_data,
    const std::vector<double>& recoil_energy_grid,
    const unsigned subshell,
    const double binding_energy,
	std::shared_ptr<const ElectroionizationSubshellElectronScatteringDistribution>&
	  electroionization_subshell_distribution )
{
  // Subshell distribution
  ElectroionizationSubshellElectronScatteringDistribution::ElectroionizationSubshellDistribution
        subshell_distribution( recoil_energy_grid.size() );

  // Create the subshell distribution
  createSubshellDistribution( raw_electroionization_data,
                              recoil_energy_grid,
                              subshell,
	                          subshell_distribution );

  electroionization_subshell_distribution.reset(
    new ElectroionizationSubshellElectronScatteringDistribution(
            subshell_distribution,
            binding_energy ) );
}

// Create the subshell recoil distribution
void ElectroionizationSubshellElectronScatteringDistributionENDLFactory::createSubshellDistribution(
	const Data::ENDLDataContainer& raw_electroionization_data,
    const std::vector<double> recoil_energy_grid,
    const unsigned subshell,
    ElectroionizationSubshellElectronScatteringDistribution::ElectroionizationSubshellDistribution&
	 subshell_distribution )
{
  for( unsigned n = 0; n < recoil_energy_grid.size(); ++n )
  {
    subshell_distribution[n].first = recoil_energy_grid[n];

    // Get the recoil energy distribution at the incoming energy
    Teuchos::Array<double> recoil_energy(
        raw_electroionization_data.getElectroionizationRecoilEnergyAtEnergy(
            subshell,
            recoil_energy_grid[n] ) );

    // Get the recoil energy pdf at the incoming energy
    Teuchos::Array<double> pdf(
        raw_electroionization_data.getElectroionizationRecoilPDFAtEnergy(
            subshell,
            recoil_energy_grid[n] ) );

    subshell_distribution[n].second.reset(
	  new const Utility::TabularDistribution<Utility::LinLin>( recoil_energy,
                                                               pdf ) );
  }
}

} // end MonteCarlo namespace

//---------------------------------------------------------------------------//
// end MonteCarlo_ElectroionizationSubshellScatteringDistributionENDLFactory.cpp
//---------------------------------------------------------------------------//

