//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_ParticleSourcePhaseSpacePoint.hpp
//! \author Alex Robinson
//! \brief  Particle source phase space point class declaration
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_PARTICLE_SOURCE_PHASE_SPACE_POINT_HPP
#define MONTE_CARLO_PARTICLE_SOURCE_PHASE_SPACE_POINT_HPP

// Std Lib Includes
#include <memory>

namespace MonteCarlo{

//! The particle source phase space point class
class ParticleSourcePhaseSpacePoint
{

public:

  //! Constructor
  ParticleSourcePhaseSpacePoint(
            const std::shared_ptr<const SpatialCoordinateConversionPolicy>&
            spatial_coord_conversion_policy,
            const std::shared_ptr<const DirectionalCoordinateConversionPolicy>&
            directional_coord_conversion_policy );

  //! Destructor
  ~ParticleSourcePhaseSpacePoint()
  { /* ... */ }

  //! Return the primary spatial coordinate of the phase space point
  double getPrimarySpatialCoordinate() const;

  //! Set the primary spatial coordinate of the phase space point
  void setPrimarySpatialCoordinate( const double primary_spatial_coord );

  //! Return the primary spatial coordinate weight
  double getPrimarySpatialCoordinateWeight() const;

  //! Set the primary spatial coordinate weight
  void setPrimarySpatialCoordinateWeight( const double weight );

  //! Return the secondary spatial coordinate of the phase space point
  double getSecondarySpatialCoordinate() const;

  //! Set the secondary spatial coordinate of the phase space point
  void setSecondarySpatialCoordinate( const double secondary_spatial_coord );

  //! Return the secondary spatial coordinate weight
  double getSecondarySpatialCoordinateWeight() const;

  //! Set the secondary spatial coordinate weight
  void setSecondarySpatialCoordinateWeight( const double weight );

  //! Return the tertiary spatial coordinate of the phase space point
  double getTertiarySpatialCoordinate() const;

  //! Set the tertiary spatial coordinate of the phase space point
  void setTertiarySpatialCoordinate( const double tertiary_spatial_coord );

  //! Return the tertiary spatial coordinate weight
  double getTertiarySpatialCoordinateWeight() const;

  //! Set the tertiary spatial coordinate weight
  void setTertiarySpatialCoordinateWeight( const double weight );

  //! Convert spatial coordinates to cartesian coordinates
  void convertSpatialCoordinatesToCartesianCoordinates(
                                               double& x_spatial_coord,
                                               double& y_spatial_coord,
                                               double& z_spatial_coord ) const;

  //! Return the weight of all spatial coordinates
  void getWeightOfSpatialCoordinates() const;

  //! Return the polar angle directional coordinate of the phase space point
  double getPolarAngleDirectionalCoordinate() const;

  //! Set the polar angle directional coordinate of the phase space point
  void setPolarAngleDirectionalCoordinate(
                                  const double polar_angle_directional_coord );

  //! Return the polar angle directional coordinate weight
  double getPolarAngleDirectionalCoordinateWeight() const;

  //! Set the polar angle directional coordinate weight
  void setPolarAngleDirectionalCoordinateWeight( const double weight );

  //! Return the azimuthal angle directional coord. of the phase space point
  double getAzimuthalAngleDirectionalCoordinate() const;

  //! Set the the azimuthal angle directional coord. of the phase space point
  void setAzimuthalAngleDirectionalCoordinate(
                              const double azimuthal_angle_directional_coord );

  //! Return the azimuthal angle directional coordinate weight
  double getAzimuthalAngleDirectionalCoordinateWeight() const;

  //! Set the azimuthal angle directional coordinate weight
  void setAzimuthalAngleDirectionalCoordinateWeight( const double weight );

  //! Convert directional coordinates to cartesian coordinates
  void convertDirectionalCoordinatesToCartesianCoordinates(
                                           double& x_directional_coord,
                                           double& y_directional_coord,
                                           double& z_directional_coord ) const;

  //! Return the weight of all directional coordinates
  double getWeightOfDirectionalCoordinates() const;

  //! Return the energy coordinate of the phase space point
  double getEnergyCoordinate() const;

  //! Set the energy coordinate of the phase space point
  void setEnergyCoordinate( const double energy_coord );

  //! Return the energy coordinate weight
  double getEnergyCoordinateWeight() const;

  //! Set the energy coordinate weight
  void setEnergyCoordinateWeight( const double weight );

  //! Return the time coordinate of the phase space point
  double getTimeCoordinate() const;

  //! Set the time coordinate of the phase space point
  void setTimeCoordinate( const double time_coord );

  //! Return the time coordinate weight
  double getTimeCoordinateWeight() const;

  //! Set the time coordinate weight
  void setTimeCoordinateWeight( const double weight );

  //! Return the weight coordinate of the phase space point
  void getWeightCoordinate() const;

  //! Set the weight coordinate of the phase space point
  void setWeightCoordinate( const double weight_coord );

  //! Return the weight of all coordinates
  void getWeightOfCoordinates() const;

private:

  // The spatial coordinate conversion policy
  std::shared_ptr<const SpatialCoordinateConversionPolicy>
  d_spatial_coord_conversion_policy;

  // The directional coordinate conversion policy
  std::shared_ptr<const DirectionalCoordinateConversionPolicy>
  d_directional_coord_conversion_policy;

  // The primary spatial coordinate
  double d_primary_spatial_coord;

  // The primary spatial coordinate weight
  double d_primary_spatial_coord_weight;

  // The secondary spatial coordinate
  double d_secondary_spatial_coord;

  // The secondary spatial coordinate weight
  double d_secondary_spatial_coord_weight;

  // The tertiary spatial coordinate
  double d_tertiary_spatial_coord;

  // The tertiary spatial coordinate weight
  double d_tertiary_spatial_coord_weight;

  // The polar angle directional coordinate
  double d_polar_angle_directional_coordinate;

  // The polar angle directional coordinate weight
  double d_polar_angle_directional_coordinate_weight;

  // The azimuthal angle directional coordinate
  double d_azimuthal_angle_directional_coordinate;

  // The azimuthal angle directional coordinate weight
  double d_azimuthal_angle_directional_coordinate_weight;

  // The energy coordinate
  double d_energy_coord;

  // The energy coordinate weight
  double d_energy_coord_weight;

  // The time coordinate
  double d_time_coord;

  // The time coordinate weight
  double d_time_coord_weight;

  // The weight coordinate
  double d_weight_coord;
};
  
} // end MonteCarlo namespace

#endif // end MONTE_CARLO_PARTICLE_SOURCE_PHASE_SPACE_POINT_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_ParticleSourcePhaseSpacePoint.hpp
//---------------------------------------------------------------------------//
