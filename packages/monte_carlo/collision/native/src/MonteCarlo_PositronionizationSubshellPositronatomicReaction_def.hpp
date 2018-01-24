//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_PositronionizationSubshellPositronatomicReaction_def.hpp
//! \author Luke Kersting
//! \brief  The positron-ionization subshell positron-atomic reaction class definition
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_POSITRONIONIZATION_SUBSHELL_POSITRONATOMIC_REACTION_DEF_HPP
#define MONTE_CARLO_POSITRONIONIZATION_SUBSHELL_POSITRONATOMIC_REACTION_DEF_HPP

// Boost Includes
#include <boost/function.hpp>
#include <boost/bind.hpp>

// FRENSIE Includes
#include "MonteCarlo_PositronionizationPositronatomicReaction.hpp"
#include "MonteCarlo_PositronatomicReactionType.hpp"
#include "Utility_ContractException.hpp"

namespace MonteCarlo{

// Basic Constructor
template<typename InterpPolicy, bool processed_cross_section>
PositronionizationSubshellPositronatomicReaction<InterpPolicy,processed_cross_section>::PositronionizationSubshellPositronatomicReaction(
    const Teuchos::ArrayRCP<const double>& incoming_energy_grid,
    const Teuchos::ArrayRCP<const double>& cross_section,
    const unsigned threshold_energy_index,
    const Data::SubshellType interaction_subshell,
    const std::shared_ptr<const ElectroionizationSubshellElectronScatteringDistribution>&
            electroionization_subshell_distribution,
    const double min_electron_energy )
  : BaseType( incoming_energy_grid,
              cross_section,
              threshold_energy_index ),
    d_interaction_subshell( interaction_subshell ),
    d_electroionization_subshell_distribution( electroionization_subshell_distribution ),
    d_reaction_type( convertSubshellEnumToPositronionizationPositronatomicReactionEnum(
                                                        interaction_subshell ) )
{
  // Make sure the interaction subshell is valid
  testPrecondition( interaction_subshell != Data::INVALID_SUBSHELL );
  testPrecondition( interaction_subshell !=Data::UNKNOWN_SUBSHELL );

  // Make sure the distribution data is valid
  testPrecondition( electroionization_subshell_distribution.use_count() > 0 );

  // Make sure the threshold energy isn't less than the binding energy
  testPrecondition( incoming_energy_grid[threshold_energy_index] >=
                    d_electroionization_subshell_distribution->getBindingEnergy() );

  // Make sure the min electron energy is valid
  testPrecondition( min_electron_energy > 0.0 );
  testPrecondition( min_electron_energy <
                    incoming_energy_grid[incoming_energy_grid.size()-1] );

  // Set the min electron energy
  this->setMinElectronEnergy( min_electron_energy );
}


// Constructor
template<typename InterpPolicy, bool processed_cross_section>
PositronionizationSubshellPositronatomicReaction<InterpPolicy,processed_cross_section>::PositronionizationSubshellPositronatomicReaction(
    const Teuchos::ArrayRCP<const double>& incoming_energy_grid,
    const Teuchos::ArrayRCP<const double>& cross_section,
    const unsigned threshold_energy_index,
    const Teuchos::RCP<const Utility::HashBasedGridSearcher>& grid_searcher,
    const Data::SubshellType interaction_subshell,
    const std::shared_ptr<const ElectroionizationSubshellElectronScatteringDistribution>&
            electroionization_subshell_distribution,
    const double min_electron_energy )
  : BaseType( incoming_energy_grid,
              cross_section,
              threshold_energy_index,
              grid_searcher ),
    d_interaction_subshell( interaction_subshell ),
    d_electroionization_subshell_distribution(
            electroionization_subshell_distribution ),
    d_reaction_type( convertSubshellEnumToPositronionizationPositronatomicReactionEnum(
            interaction_subshell ) )
{
  // Make sure the interaction subshell is valid
  testPrecondition( interaction_subshell != Data::INVALID_SUBSHELL );
  testPrecondition( interaction_subshell !=Data::UNKNOWN_SUBSHELL );

  // Make sure the distribution data is valid
  testPrecondition( electroionization_subshell_distribution.use_count() > 0 );

  // Make sure the threshold energy isn't less than the binding energy
  testPrecondition( incoming_energy_grid[threshold_energy_index] >=
                    d_electroionization_subshell_distribution->getBindingEnergy() );

  // Make sure the min electron energy is valid
  testPrecondition( min_electron_energy > 0.0 );
  testPrecondition( min_electron_energy <
                    incoming_energy_grid[incoming_energy_grid.size()-1] );

  // Set the min electron energy
  this->setMinElectronEnergy( min_electron_energy );
}

// Return the differential cross section
template<typename InterpPolicy, bool processed_cross_section>
double PositronionizationSubshellPositronatomicReaction<InterpPolicy,processed_cross_section>::getDifferentialCrossSection(
    const double incoming_energy,
    const double outgoing_energy ) const
{
  // Make sure the energies are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( outgoing_energy >= 0.0 );
  testPrecondition( outgoing_energy <= incoming_energy );

  if ( !this->isEnergyWithinEnergyGrid( incoming_energy ) )
    return 0.0;

  // Evaluate the forward cross section at the incoming energy
  double forward_cs = this->getCrossSection( incoming_energy );

  // Sample the pdf using the energy of the knock-on electron
  double pdf = d_electroionization_subshell_distribution->evaluatePDF(
          incoming_energy,
          outgoing_energy );

  return forward_cs*pdf;
}

// Simulate the reaction
template<typename InterpPolicy, bool processed_cross_section>
void PositronionizationSubshellPositronatomicReaction<InterpPolicy,processed_cross_section>::react(
     PositronState& positron,
     ParticleBank& bank,
     Data::SubshellType& shell_of_interaction ) const
{
  // Make sure the positron energy isn't less than the binding energy
  testPrecondition( positron.getEnergy() >=
                    d_electroionization_subshell_distribution->getBindingEnergy() );

  d_electroionization_subshell_distribution->scatterPositron(
                                               positron,
                                               bank,
                                               shell_of_interaction);

  positron.incrementCollisionNumber();

  shell_of_interaction = d_interaction_subshell;

  // Annihilate the positron if it is below cutoff
  if ( positron.getEnergy() < this->getMinElectronEnergy() )
  {
    this->annihilatePositron( positron, bank );
  }
}

// Return the reaction type
template<typename InterpPolicy, bool processed_cross_section>
PositronatomicReactionType PositronionizationSubshellPositronatomicReaction<InterpPolicy,processed_cross_section>::getReactionType() const
{
  return d_reaction_type;
}

// Get the interaction subshell (non-standard interface)
template<typename InterpPolicy, bool processed_cross_section>
unsigned PositronionizationSubshellPositronatomicReaction<InterpPolicy,processed_cross_section>::getSubshell() const
{
  return d_interaction_subshell;
}

} // end MonteCarlo namespace

#endif // end MONTE_CARLO_POSITRONIONIZATION_SUBSHELL_POSITRONATOMIC_REACTION_DEF_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_PositronionizationSubshellPositronatomicReaction_def.hpp
//---------------------------------------------------------------------------//
