//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_SubshellIncoherentPhotonScatteringDistribution.hpp
//! \author Alex Robinson
//! \brief  The subshell incoherent photon scattering distribution decl.
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_SUBSHELL_INCOHERENT_PHOTON_SCATTERING_DISTRIBUTION_HPP
#define MONTE_CARLO_SUBSHELL_INCOHERENT_PHOTON_SCATTEIRNG_DISTRIBUTION_HPP

// Boost Includes
#include <boost/function.hpp>

// Trilinos Includes
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>

// FRENSIE Includes
#include "MonteCarlo_PhotonScatteringDistribution.hpp"
#include "MonteCarlo_SubshellType.hpp"
#include "Utility_TabularOneDDistribution.hpp"

namespace MonteCarlo{

//! The subshell incoherent photon scattering distribution class
class SubshellIncoherentPhotonScatteringDistribution : public PhotonScatteringDistribution
{

public:

  //! Constructor without Doppler broadening
  SubshellIncoherentPhotonScatteringDistribution(
      const SubshellType interaction_subshell,
      const double binding_energy,
      const Teuchos::RCP<const Utility::OneDDistribution>& occupation_number );

  //! Constructor for Doppler broadening
  SubshellIncoherentPhotonScatteringDistribution(
	const SubshellType interaction_subshell,
	const double binding_energy,
	const Teuchos::RCP<const Utility::OneDDistribution>& occupation_number,
        const Teuchos::RCP<const Utility::TabularOneDDistribution>& 
	compton_profile );

  //! Destructor
  ~SubshellIncoherentPhotonScatteringDistribution()
  { /* ... */ }

  //! Return the subshell
  SubshellType getSubshell() const;

  //! Return the binding energy
  double getBindingEnergy() const;

  //! Evaluate the distribution
  double evaluate( const double incoming_energy,
		   const double scattering_angle_cosine ) const;

  //! Evaluate the PDF
  double evaluatePDF( const double incoming_energy,
		      const double scattering_angle_cosine ) const;
  
  //! Evaluate the integrated cross section (cm^2)
  double evaluateIntegratedCrossSection( const double incoming_energy,
					 const double precision ) const;

  //! Sample an outgoing energy and direction from the distribution
  void sample( const double incoming_energy,
	       double& outgoing_energy,
	       double& scattering_angle_cosine,
	       SubshellType& shell_of_interaction ) const;

  //! Sample an outgoing energy and direction and record the number of trials
  void sampleAndRecordTrials( const double incoming_energy,
			      double& outgoing_energy,
			      double& scattering_angle_cosine,
			      SubshellType& shell_of_interaction,
			      unsigned& trials ) const;

  //! Randomly scatter the photon
  void scatterPhoton( PhotonState& photon,
		      ParticleBank& bank,
		      SubshellType& shell_of_interaction ) const;

private:

  // Ignore Doppler broadening
  double returnComptonLine( 
			 const double initial_energy,
			 const double compton_line_energy,
			 const double scattering_angle_cosine,
			 const double max_electron_momentum_projection ) const;
			    
  // Doppler broaden a compton line
  double dopplerBroadenComptonLine(
			 const double initial_energy,
			 const double compton_line_energy,
			 const double scattering_angle_cosine,
			 const double max_electron_momentum_projection ) const;

  // The interaction subshell
  SubshellType d_subshell;

  // The subshell binding energy
  double d_binding_energy;

  // The occupation number distribution
  Teuchos::RCP<const Utility::OneDDistribution> d_occupation_number;

  // The compton profile distribution
  Teuchos::RCP<const Utility::TabularOneDDistribution> d_compton_profile;

  // The Doppler broadening function pointer
  boost::function<double (double, double, double, double)>
  d_doppler_broadening_func;
};

} // end MonteCarlo namespace

#endif // end MONTE_CARLO_SUBSHELL_INCOHERENT_PHOTON_SCATTERING_DISTRIBUTION_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_SubshellIncoherentPhotonScatteringDistribution.hpp
//---------------------------------------------------------------------------//
