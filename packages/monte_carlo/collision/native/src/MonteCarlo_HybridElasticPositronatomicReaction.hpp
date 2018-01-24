//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_HybridElasticPositronatomicReaction.hpp
//! \author Luke Kersting
//! \brief  The hybrid scattering elastic positron-atomic reaction class decl.
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_HYBRID_ELASTIC_POSITRONATOMIC_REACTION_HPP
#define MONTE_CARLO_HYBRID_ELASTIC_POSITRONATOMIC_REACTION_HPP

// Trilinos Includes
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "MonteCarlo_PositronatomicReaction.hpp"
#include "MonteCarlo_StandardGenericAtomicReaction.hpp"
#include "MonteCarlo_HybridElasticElectronScatteringDistribution.hpp"

namespace MonteCarlo{

//! The hybrid elastic positron-atomic reaction class
template<typename InterpPolicy,
         bool processed_cross_section = false>
class HybridElasticPositronatomicReaction : public StandardGenericAtomicReaction<PositronatomicReaction,InterpPolicy,processed_cross_section>
{

private:

  // Typedef for the base class type
typedef StandardGenericAtomicReaction<PositronatomicReaction,InterpPolicy,processed_cross_section>
    BaseType;

public:

  //! Basic Constructor
  HybridElasticPositronatomicReaction(
    const Teuchos::ArrayRCP<const double>& incoming_energy_grid,
    const Teuchos::ArrayRCP<const double>& cross_section,
    const unsigned threshold_energy_index,
    const double cutoff_angle_cosine,
    const std::shared_ptr<const HybridElasticElectronScatteringDistribution>&
            hybrid_distribution );

  //! Constructor
  HybridElasticPositronatomicReaction(
    const Teuchos::ArrayRCP<const double>& incoming_energy_grid,
    const Teuchos::ArrayRCP<const double>& cross_section,
    const unsigned threshold_energy_index,
    const Teuchos::RCP<const Utility::HashBasedGridSearcher>& grid_searcher,
    const double cutoff_angle_cosine,
    const std::shared_ptr<const HybridElasticElectronScatteringDistribution>&
            hybrid_distribution );


  //! Destructor
  ~HybridElasticPositronatomicReaction()
  { /* ... */ }

  //! Return the number of electrons emitted from the rxn at the given energy
  unsigned getNumberOfEmittedElectrons( const double energy ) const;

  //! Return the number of photons emitted from the rxn at the given energy
  unsigned getNumberOfEmittedPhotons( const double energy ) const;

  //! Return the reaction type
  PositronatomicReactionType getReactionType() const;

  //! Return the differential cross section
  double getDifferentialCrossSection( const double incoming_energy,
                                      const double scattering_angle_cosine ) const;

  //! Simulate the reaction
  void react( PositronState& positron,
              ParticleBank& bank,
              Data::SubshellType& shell_of_interaction ) const;

private:

  // The hybrid elastic scattering distribution
  std::shared_ptr<const HybridElasticElectronScatteringDistribution>
    d_hybrid_distribution;
};

} // end MonteCarlo namespace

//---------------------------------------------------------------------------//
// Template Includes
//---------------------------------------------------------------------------//

#include "MonteCarlo_HybridElasticPositronatomicReaction_def.hpp"

//---------------------------------------------------------------------------//

#endif // end MONTE_CARLO_HYBRID_ELASTIC_POSITRONATOMIC_REACTION_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_HybridElasticPositronatomicReaction.hpp
//---------------------------------------------------------------------------//
