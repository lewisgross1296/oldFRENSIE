//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_TetMeshTrackLengthFluxEstimator.hpp
//! \author Alex Robinson, Eli Moll
//! \brief  Tet mesh flux estimator class declaration.
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_TET_MESH_TRACK_LENGTH_FLUX_ESTIMATOR_HPP
#define MONTE_CARLO_TET_MESH_TRACK_LENGTH_FLUX_ESTIMATOR_HPP

// Std Lib Includes
#include <string>

// Boost Includes
#include <boost/mpl/vector.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

// Moab Includes
#include <moab/Interface.hpp>
#include <moab/AdaptiveKDTree.hpp>
#include <moab/OrientedBoxTreeTool.hpp>
#include <moab/Matrix3.hpp>

// Trilinos Includes
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "MonteCarlo_StandardEntityEstimator.hpp"
#include "MonteCarlo_ParticleSubtrackEndingGlobalEventObserver.hpp"
#include "MonteCarlo_EstimatorContributionMultiplierPolicy.hpp"
#include "Geometry_ModuleTraits.hpp"
#include "MonteCarlo_ParticleState.hpp"

namespace MonteCarlo{

//! The tet-mesh track length flux estimator class
template<typename ContributionMutliplierPolicy = WeightMultiplier>
class TetMeshTrackLengthFluxEstimator : public StandardEntityEstimator<moab::EntityHandle>,
  public ParticleSubtrackEndingGlobalEventObserver
{

private:

  // Typedef for triangle intersection pairs
  typedef Utility::Pair<double,moab::EntityHandle> IntersectionData;

public:

  //! Typedef for the cell id type
  typedef Geometry::ModuleTraits::InternalCellHandle cellIdType;

  //! Typedef for event tags used for quick dispatcher registering
  typedef boost::mpl::vector<ParticleSubtrackEndingGlobalEventObserver::EventTag>
  EventTags;

  //! Constructor
  TetMeshTrackLengthFluxEstimator(
		     const Estimator::idType id,
		     const double multiplier,
		     const std::string input_mesh_file_name,
		     const std::string output_mesh_file_name = "tetmesh.h5m" );

  //! Destructor
  ~TetMeshTrackLengthFluxEstimator()
  { /* ... */ }

  //! Set the response functions
  void setResponseFunctions(
  const Teuchos::Array<Teuchos::RCP<ResponseFunction> >& response_functions );

  //! Set the particle types that can contribute to the estimator
  void setParticleTypes( const Teuchos::Array<ParticleType>& particle_types );

  //! Add current history estimator contribution
  void updateFromGlobalParticleSubtrackEndingEvent(
						 const ParticleState& particle,
						 const double start_point[3],
						 const double end_point[3] );

  //! Export the estimator data
  void exportData( EstimatorHDF5FileHandler& hdf5_file,
		   const bool process_data ) const;

  //! Print the estimator data
  void print( std::ostream& os ) const;

  //! Get all tet elements
  const moab::Range getAllTetElements() const;

  //! Determine which tet the point is in
  moab::EntityHandle whichTetIsPointIn( const double point[3] );

private:

  // Compare intersection data
  static bool compareIntersections( const IntersectionData& a,
				    const IntersectionData& b );

  // Assign bin boundaries to an estimator dimension
  void assignBinBoundaries(
	const Teuchos::RCP<EstimatorDimensionDiscretization>& bin_boundaries );

  // The moab instance that stores all mesh data
  Teuchos::RCP<moab::Interface> d_moab_interface;

  // The tet meshset
  moab::EntityHandle d_tet_meshset;

  // The kd-tree for finding point in tet
  Teuchos::RCP<moab::AdaptiveKDTree> d_kd_tree;
  
  // The root of the kd-tree
  moab::EntityHandle d_kd_tree_root;
  
  // The oriented box tree for finding intersections with the mesh
  Teuchos::RCP<moab::OrientedBoxTreeTool> d_obb_tree;

  // The root of the obb-tree
  moab::EntityHandle d_obb_tree_root;

  // The last cell that was visited
  Geometry::ModuleTraits::InternalCellHandle d_last_visited_cell;

  // The map of tet ids and barycentric coordinate transform matrices
  boost::unordered_map<moab::EntityHandle,moab::Matrix3> 
  d_tet_barycentric_transform_matrices;
  
  // The map of tet ids and reference vertices
  boost::unordered_map<moab::EntityHandle, moab::CartVect>
  d_tet_reference_vertices;
  
  // The output mesh file name
  std::string d_output_mesh_name;
};
  
} // end MonteCarlo namespace

//---------------------------------------------------------------------------//
// Template Includes
//---------------------------------------------------------------------------//

#include "MonteCarlo_TetMeshTrackLengthFluxEstimator_def.hpp"

//---------------------------------------------------------------------------//

#endif // end MONTE_CARLO_TET_MESH_TRACK_LENGTH_FLUX_ESTIMATOR_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_TetMeshTrackLengthFluxEstimator.hpp
//---------------------------------------------------------------------------//
