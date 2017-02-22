//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_Estimator.cpp
//! \author Alex Robinson
//! \brief  Estimator base class definition
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <limits>

// FRENSIE Includes
#include "MonteCarlo_Estimator.hpp"
#include "MonteCarlo_EstimatorHDF5FileHandler.hpp"

namespace MonteCarlo{

// Initialize the tolerance
double Estimator::tol = 1e-8;

// Constructor
  Estimator::Estimator( const ParticleHistoryObserver::idType id,
			const double multiplier )
    : ParticleHistoryObserver( id ),
      d_multiplier( multiplier ),
      d_has_uncommitted_history_contribution( 1, false ),
      d_response_functions( 1 )
{
  // Make sure the multiplier is valid
  testPrecondition( multiplier > 0.0 );

  // Set the response function
  d_response_functions[0] = ResponseFunction::default_response_function;
}

// Check if the estimator has uncommitted history contributions
bool Estimator::hasUncommittedHistoryContribution(
					       const unsigned thread_id ) const
{
  // Make sure the thread is is valid
  testPrecondition( thread_id < d_has_uncommitted_history_contribution.size());

  return d_has_uncommitted_history_contribution[thread_id];
}

// Enable support for multiple threads
void Estimator::enableThreadSupport( const unsigned num_threads )
{
  // Make sure only the master thread calls this function
  testPrecondition( Utility::GlobalOpenMPSession::getThreadId() == 0 );

  d_has_uncommitted_history_contribution.resize( num_threads, false );
}

// Export the estimator data
void Estimator::exportData(
                    const std::shared_ptr<Utility::HDF5FileHandler>& hdf5_file,
                    const bool process_data ) const
{
  // Make sure only the master thread calls this function
  testPrecondition( Utility::GlobalOpenMPSession::getThreadId() == 0 );
  // Make sure this estimator has not been exported yet
  remember( EstimatorHDF5FileHandler test_estimator_hdf5_file( hdf5_file ) );
  testPrecondition(
               !test_estimator_hdf5_file.doesEstimatorExist( this->getId() ) );

  // Open the estimator hdf5 file
  EstimatorHDF5FileHandler estimator_hdf5_file( hdf5_file );
  
  // Export the estimator multiplier
  estimator_hdf5_file.setEstimatorMultiplier( this->getId(), d_multiplier );

  // Export the response function ordering
  {
    Teuchos::Array<size_t> response_function_ordering(
						 d_response_functions.size() );
    for( size_t i = 0; i < d_response_functions.size(); ++i )
      response_function_ordering[i] = d_response_functions[i]->getId();

    estimator_hdf5_file.setEstimatorResponseFunctionOrdering(
						  this->getId(),
						  response_function_ordering );
  }

  // Export the dimension discretization
  d_dimension_discretization.export( this->getId(), estimator_hdf5_file );
}

// Assign response function to the estimator
/*! \details Override this method in a derived class if response function
 * properties need to be checked before the assignment takes place.
 */
void Estimator::assignResponseFunction(
                             const ResponseFunctionPointer& response_function )
{
  // Make sure only the master thread calls this function
  testPrecondition( Utility::GlobalOpenMPSession::getThreadId() == 0 );
  // Make sure that the response function pointer is valid
  testPrecondition( response_function.get() );
  
  d_response_functions.push_back( response_function );
}

// Assign the particle type to the estimator
/*! \details Override this method in a derived class if response function
 * properties need to be checked before assignment takes place.
 */
void Estimator::assignParticleType( const ParticleType particle_type )
{
  // Make sure only the master thread calls this function
  testPrecondition( Utility::GlobalOpenMPSession::getThreadId() == 0 );

  d_particle_types.insert( particle_type );
}

// Assign bins to an estimator dimension
void Estimator::assignBins( const DimensionDiscretizationPoint& discretization )
{
  // Make sure only the master thread calls this function
  testPrecondition( Utility::GlobalOpenMPSession::getThreadId() == 0 );

  d_phase_space_discretization.assignDiscretizationToDimension( discretization );
}

// Set the has uncommited history contribution flag
/*! \details This should be called whenever the current history contributes
 * to the estimator.
 */
void Estimator::setHasUncommittedHistoryContribution(
						     const unsigned thread_id )
{
  // Make sure the thread is is valid
  testPrecondition( thread_id < d_has_uncommitted_history_contribution.size());

  d_has_uncommitted_history_contribution[thread_id] = true;
}

// Unset the has uncommited history contribution flag
/*! \details This should be called when the current history contribution is
 * committed to the estimator
 */
void Estimator::unsetHasUncommittedHistoryContribution(
						     const unsigned thread_id )
{
  // Make sure the thread is is valid
  testPrecondition( thread_id < d_has_uncommitted_history_contribution.size());

  d_has_uncommitted_history_contribution[thread_id] = false;
}

// Return the bin name
std::string Estimator::getBinName( const size_t bin_index ) const
{
  // Make sure the bin index is valid
  testPrecondition( bin_index <
		    getNumberOfBins()*getNumberOfResponseFunctions() );

  return d_phase_space_discretization.getBinName( bin_index ) +
    this->getResponseFunctionName(
                           this->calculateResponseFunctionIndex( bin_index ) );
}

// Print the estimator response function names
void Estimator::printEstimatorResponseFunctionNames( std::ostream& os ) const
{
  os << "Response Functions: " << std::endl;

  for( size_t i = 0u; i < d_response_functions.size(); ++i )
    os << "  " << i+1 << ".) " << getResponseFunctionName( i ) << std::endl;
}

// Print the estimator bins
void Estimator::printEstimatorBins( std::ostream& os ) const
{
  for( size_t i = 0u; i < d_dimension_ordering.size(); ++i )
  {
    d_dimension_bins_map.find(d_dimension_ordering[i])->second->print( os );
    os << std::endl;
  }
}

// Print the estimator data stored in an array
/*! \details The number of elements in the array should be equal to the
 * the number of estimator bins times the number of response functions.
 */
void Estimator::printEstimatorBinData(
			std::ostream& os,
		        const TwoEstimatorMomentsArray& estimator_moments_data,
			const double norm_constant ) const
{
  // Make sure that the estimator moment array is valid
  testPrecondition( estimator_moments_data.size() ==
		    getNumberOfBins()*getNumberOfResponseFunctions() );

  // Use this array to determine when bin indices should be printed
  Teuchos::Array<size_t> previous_bin_indices(
					d_dimension_ordering.size(),
					std::numeric_limits<size_t>::max() );

  for( size_t r = 0u; r < getNumberOfResponseFunctions(); ++r )
  {
    os << "Response Function: " << getResponseFunctionName( r ) << std::endl;

    for( size_t i = 0u; i < getNumberOfBins(); ++i )
    {
      for( int d = d_dimension_ordering.size()-1; d >= 0; --d )
      {
	const size_t& dimension_index_step_size =
	 d_dimension_index_step_size_map.find(d_dimension_ordering[d])->second;

	// Calculate the bin index for the dimension
	size_t bin_index = (i/dimension_index_step_size)%
	  (getNumberOfBins( d_dimension_ordering[d] ));

	// Print the bin boundaries if the dimension index has changed
	if( bin_index != previous_bin_indices[d] )
	{
	  previous_bin_indices[d] = bin_index;

	  // Add proper spacing before printing
	  for( size_t s = 0u; s < d_dimension_ordering.size()-d; ++s )
	    os << " ";

	  d_dimension_bins_map.find(d_dimension_ordering[d])->second->printBoundariesOfBin( os, bin_index );

	  // Print a new line character for all but the first dimension
	  if( d != 0 )
	    os << std::endl;
	}
      }

      // Calculate the bin index for the response function
      size_t bin_index = i + r*getNumberOfBins();


      // Calculate the estimator bin data
      double estimator_bin_value;
      double estimator_bin_rel_err;

      processMoments( estimator_moments_data[bin_index],
		      norm_constant,
		      estimator_bin_value,
		      estimator_bin_rel_err );

      // Print the estimator bin data
      os << " " << estimator_bin_value << " "
	 << estimator_bin_rel_err << std::endl;
    }
  }
}

// Print the total estimator data stored in an array
/*! \details The array should have one element for each response function
 * assigned to the estimator.
 */
void Estimator::printEstimatorTotalData(
		 std::ostream& os,
		 const FourEstimatorMomentsArray& total_estimator_moments_data,
		 const double norm_constant ) const
{
  // Make sure that the total estimator moments data is valid
  testPrecondition( total_estimator_moments_data.size() ==
		    getNumberOfResponseFunctions() );

   for( size_t i = 0; i < d_response_functions.size(); ++i )
  {
    os << "Response Function: " << getResponseFunctionName( i ) << std::endl;

    double estimator_value;
    double estimator_rel_err;
    double estimator_vov;
    double estimator_fom;

    processMoments( total_estimator_moments_data[i],
		    norm_constant,
		    estimator_value,
		    estimator_rel_err,
		    estimator_vov,
		    estimator_fom );

    os << estimator_value << " "
       << estimator_rel_err << " "
       << estimator_vov << " "
       << estimator_fom << std::endl;
  }
}

// Calculate the response function index given a bin index
size_t Estimator::calculateResponseFunctionIndex(
					       const size_t bin_index ) const
{
  // Make sure the bin index is valid
  testPrecondition( bin_index <
		    getNumberOfBins()*getNumberOfResponseFunctions() );
  testPrecondition( bin_index < std::numeric_limits<size_t>::max() );

  return bin_index/getNumberOfBins();
}

// Convert first and second moments to mean and relative error
void Estimator::processMoments( const Utility::Pair<double,double>& moments,
				const double norm_constant,
				double& mean,
				double& relative_error ) const
{
  // Make sure the moments are valid
  testPrecondition( !ST::isnaninf( moments.first ) );
  testPrecondition( !ST::isnaninf( moments.second ) );
  // Make sure the norm contant is valid
  testPrecondition( !ST::isnaninf( norm_constant ) );
  testPrecondition( norm_constant > 0.0 );

  mean = calculateMean( moments.first )*d_multiplier/norm_constant;

  relative_error = calculateRelativeError( moments.first, moments.second );
}

// Convert first, second, third, fourth moments to mean, rel. er., vov, fom
void Estimator::processMoments(
		     const Utility::Quad<double,double,double,double>& moments,
		     const double norm_constant,
		     double& mean,
		     double& relative_error,
		     double& variance_of_variance,
		     double& figure_of_merit ) const
{
  // Make sure the moments are valid
  testPrecondition( !ST::isnaninf( moments.first ) );
  testPrecondition( !ST::isnaninf( moments.second ) );
  testPrecondition( !ST::isnaninf( moments.third ) );
  testPrecondition( !ST::isnaninf( moments.fourth ) );
  // Make sure the norm contant is valid
  testPrecondition( !ST::isnaninf( norm_constant ) );
  testPrecondition( norm_constant > 0.0 );

  mean = calculateMean( moments.first )*d_multiplier/norm_constant;

  relative_error = calculateRelativeError( moments.first, moments.second );

  variance_of_variance = calculateVOV( moments.first,
				       moments.second,
				       moments.third,
				       moments.fourth );

  figure_of_merit = calculateFOM( relative_error );
}

// Calculate the mean of a set of contributions
double Estimator::calculateMean( const double first_moment_contributions) const
{
  // Make sure the first moment contributions are valid
  testPrecondition( !ST::isnaninf( first_moment_contributions ) );
  testPrecondition( first_moment_contributions >= 0.0 );

  double mean = first_moment_contributions/this->getNumberOfHistories();

  // Make sure the mean is valid
  testPostcondition( !ST::isnaninf( mean ) );
  testPostcondition( mean >= 0.0 );

  return mean;
}

// Calculate the relative error of a set of contributions
double Estimator::calculateRelativeError(
			       const double first_moment_contributions,
			       const double second_moment_contributions ) const
{
  // Make sure that the number of histories is valid
  testPrecondition( this->getNumberOfHistories() > 0ull );
  // Make sure the first moment contributions are valid
  testPrecondition( !ST::isnaninf( first_moment_contributions ) );
  testPrecondition( first_moment_contributions >= 0.0 );
  // Make sure the second moment contributions are valid
  testPrecondition( !ST::isnaninf( second_moment_contributions ) );
  testPrecondition( second_moment_contributions >= 0.0 );

  double relative_error;

  if( first_moment_contributions > 0.0 )
  {
    double argument = second_moment_contributions/
      (first_moment_contributions*first_moment_contributions) -
      1.0/this->getNumberOfHistories();

    // Check for roundoff error resulting in a very small negative number
    if( argument < 0.0 && argument > -Estimator::tol )
      relative_error = 0.0;
    else
      relative_error = ST::squareroot( argument );
  }
  else
    relative_error = 0.0;

  // Make sure the relative error is valid
  testPostcondition( !ST::isnaninf( relative_error ) );
  testPostcondition( relative_error >= 0.0 );

  return relative_error;
}

// Calculate the variance of the variance (VOV) of a set of contributions
/* \details VOV = S^2(S^2_xbar)/S^4_xbar
 */
double Estimator::calculateVOV( const double first_moment_contributions,
				const double second_moment_contributions,
				const double third_moment_contributions,
				const double fourth_moment_contributions) const
{
  // Make sure that the number of histories is valid
  testPrecondition( this->getNumberOfHistories() > 0ull );
  // Make sure the first moment contributions are valid
  testPrecondition( !ST::isnaninf( first_moment_contributions ) );
  testPrecondition( first_moment_contributions >= 0.0 );
  // Make sure the second moment contributions are valid
  testPrecondition( !ST::isnaninf( second_moment_contributions ) );
  testPrecondition( second_moment_contributions >= 0.0 );
  // Make sure the third moment contributions are valid
  testPrecondition( !ST::isnaninf( third_moment_contributions ) );
  testPrecondition( third_moment_contributions >= 0.0 );
  // Make sure the fourth moment contributions are valid
  testPrecondition( !ST::isnaninf( fourth_moment_contributions ) );
  testPrecondition( fourth_moment_contributions >= 0.0 );

  double first_moment_contributions_squared =
    first_moment_contributions*first_moment_contributions;

  double num_histories_squared = this->getNumberOfHistories()*
    this->getNumberOfHistories();

  double num_histories_cubed = num_histories_squared*this->getNumberOfHistories();

  double vov_numerator = fourth_moment_contributions -
    4*first_moment_contributions*third_moment_contributions/
    this->getNumberOfHistories() +
    8*second_moment_contributions*first_moment_contributions_squared/
    num_histories_squared -
    4*first_moment_contributions_squared*first_moment_contributions_squared/
    num_histories_cubed -
    second_moment_contributions*second_moment_contributions/
    this->getNumberOfHistories();

  double vov_denominator =
    (second_moment_contributions - first_moment_contributions_squared/
     this->getNumberOfHistories());
  vov_denominator *= vov_denominator;

  double vov;

  // Check for roundoff errors resulting in a very small negative number
  if( vov_denominator <= 0.0 && vov_denominator > -Estimator::tol )
    vov = 0.0;
  else if( vov_numerator <= 0.0 && vov_denominator > -Estimator::tol )
    vov = 0.0;
  else
    vov = vov_numerator/vov_denominator;

  // Make sure the variance of the variance is valid
  testPostcondition( !ST::isnaninf( vov ) );
  testPostcondition( vov >= 0.0 );

  return vov;
}

// Calculate the figure of merit (FOM) of a set of contributions
double Estimator::calculateFOM( const double relative_error ) const
{
  // Make sure the problem time is valid
  testPrecondition( this->getElapsedTime() > 0.0 );

  double problem_time = this->getElapsedTime();

  double fom;

  if( relative_error > 0.0 )
    fom = 1.0/(relative_error*relative_error*problem_time);
  else
    fom = 0.0;

  // Make sure the figure of merit is valid
  testPostcondition( !ST::isnaninf( fom ) );
  testPostcondition( fom >= 0.0 );

  return fom;
}

} // end MonteCarlo namespace

//---------------------------------------------------------------------------//
// end MonteCarlo_Estimator.cpp
//---------------------------------------------------------------------------//
