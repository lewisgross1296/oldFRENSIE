//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_CoupledElasticPositronatomicReaction_def.hpp
//! \author Luke Kersting
//! \brief  The coupled scattering elastic positron-atomic reaction class def.
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_COUPLED_ELASTIC_POSITRONATOMIC_REACTION_DEF_HPP
#define MONTE_CARLO_COUPLED_ELASTIC_POSITRONATOMIC_REACTION_DEF_HPP

// FRENSIE Includes
#include "Utility_ContractException.hpp"

namespace MonteCarlo{

// Basic Constructor
template<typename InterpPolicy, bool processed_cross_section>
CoupledElasticPositronatomicReaction<InterpPolicy,processed_cross_section>::CoupledElasticPositronatomicReaction(
      const Teuchos::ArrayRCP<const double>& incoming_energy_grid,
      const Teuchos::ArrayRCP<const double>& cross_section,
      const unsigned threshold_energy_index,
      const std::shared_ptr<const CoupledElasticElectronScatteringDistribution>&
            scattering_distribution )
  : BaseType( incoming_energy_grid,
              cross_section,
              threshold_energy_index ),
    d_scattering_distribution( scattering_distribution )
{
  // Make sure scattering distribution is valid
  testPrecondition( scattering_distribution.use_count() > 0 );
}

// Constructor
template<typename InterpPolicy, bool processed_cross_section>
CoupledElasticPositronatomicReaction<InterpPolicy,processed_cross_section>::CoupledElasticPositronatomicReaction(
      const Teuchos::ArrayRCP<const double>& incoming_energy_grid,
      const Teuchos::ArrayRCP<const double>& cross_section,
      const unsigned threshold_energy_index,
      const Teuchos::RCP<const Utility::HashBasedGridSearcher>& grid_searcher,
      const std::shared_ptr<const CoupledElasticElectronScatteringDistribution>&
            scattering_distribution )
  : BaseType( incoming_energy_grid,
              cross_section,
              threshold_energy_index,
              grid_searcher ),
    d_scattering_distribution( scattering_distribution )
{
  // Make sure scattering distribution is valid
  testPrecondition( scattering_distribution.use_count() > 0 );
}

// Return the number of photons emitted from the rxn at the given energy
/*! \details This does not include photons from atomic relaxation.
 */
template<typename InterpPolicy, bool processed_cross_section>
unsigned CoupledElasticPositronatomicReaction<InterpPolicy,processed_cross_section>::getNumberOfEmittedPhotons( const double energy ) const
{
  return 0u;
}

// Return the number of electrons emitted from the rxn at the given energy
template<typename InterpPolicy, bool processed_cross_section>
unsigned CoupledElasticPositronatomicReaction<InterpPolicy,processed_cross_section>::getNumberOfEmittedElectrons( const double energy ) const
{
  return 0u;
}

// Return the reaction type
template<typename InterpPolicy, bool processed_cross_section>
PositronatomicReactionType CoupledElasticPositronatomicReaction<InterpPolicy,processed_cross_section>::getReactionType() const
{
  return COUPLED_ELASTIC_POSITRONATOMIC_REACTION;
}

// Return the differential cross section
template<typename InterpPolicy, bool processed_cross_section>
double CoupledElasticPositronatomicReaction<InterpPolicy,processed_cross_section>::getDifferentialCrossSection(
            const double incoming_energy,
            const double scattering_angle_cosine ) const
{
  // Get the PDF
  double pdf =
    d_scattering_distribution->evaluatePDF( incoming_energy,
                                            scattering_angle_cosine );

  // Get the cross section
  double cross_section = this->getCrossSection( incoming_energy );

  return pdf*cross_section;
}

// Simulate the reaction
template<typename InterpPolicy, bool processed_cross_section>
void CoupledElasticPositronatomicReaction<InterpPolicy,processed_cross_section>::react(
                     PositronState& positron,
                     ParticleBank& bank,
                     Data::SubshellType& shell_of_interaction ) const
{
  d_scattering_distribution->scatterPositron( positron,
                                              bank,
                                              shell_of_interaction);

  positron.incrementCollisionNumber();

  // The shell of interaction is currently ignored
  shell_of_interaction =Data::UNKNOWN_SUBSHELL;
}

} // end MonteCarlo namespace

#endif // end MONTE_CARLO_COUPLED_ELASTIC_POSITRONATOMIC_REACTION_DEF_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_CoupledElasticPositronatomicReaction_def.hpp
//---------------------------------------------------------------------------//
