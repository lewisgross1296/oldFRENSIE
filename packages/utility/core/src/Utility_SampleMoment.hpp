//---------------------------------------------------------------------------//
//!
//! \file   Utility_SampleMoment.hpp
//! \author Alex Robinson
//! \brief  The sample moment class declaration.
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_SAMPLE_MOMENT_HPP
#define UTILITY_SAMPLE_MOMENT_HPP

// Std Lib Includes
#include <type_traits>

// FRENSIE Includes
#include "Utility_QuantityTraits.hpp"

namespace Utility{

/*! The sample moment class
 *
 * To avoid additional memory overhead, the number of samples
 * will not be recorded by this class. When one attempts to retrieve the 
 * raw moment the number of samples must be provided.
 */
template<size_t N, typename T, typename Enabled = void>
class SampleMoment
{ /* ... */ };

//! Partial specialization of the sample moment class for floating point types
template<size_t N, typename T>
class SampleMoment<N,T,typename std::enable_if<(N>0) && std::is_floating_point<typename QuantityTraits<T>::RawType>::value>::type>
{

public:
  
  //! The moment value type
  typedef typename QuantityTraits<T>::template GetQuantityToPowerType<N>::type
  ValueType;

  //! Process a raw score
  static ValueType processRawScore( const T& raw_score );
  
  //! Default constructor
  SampleMoment();
  
  //! Constructor
  SampleMoment( const ValueType& starting_score );

  //! Copy constructor
  SampleMoment( const SampleMoment& other_moment );

  //! Assignment operator
  SampleMoment& operator=( const SampleMoment& other_moment );

  //! Destructor
  ~SampleMoment()
  { /* ... */ }

  //! Return the current score
  const ValueType& getCurrentScore() const;

  //! Return the value of the moment given the number of samples
  ValueType get( const size_t number_of_samples ) const;

  //! Add a raw score
  void addRawScore( const T& raw_score );

  //! Add a processed score (raw_score^N)
  void addProcessedScore( const ValueType& processed_score );

private:

  // The raw scores
  ValueType d_current_score;
};

//! Calculate the mean
template<typename T>
typename SampleMoment<1,T>::ValueType calculateMean(
                                         const SampleMoment<1,T>& first_moment,
                                         const size_t number_of_samples );
  
//! Calculate the variance of the population
template<typename T>
typename SampleMoment<2,T>::ValueType calculateVariance(
                                        const SampleMoment<1,T>& first_moment,
                                        const SampleMoment<2,T>& second_moment,
                                        const size_t number_of_samples );

//! Calculate the variance of the mean
template<typename T>
typename SampleMoment<2,T>::ValueType calculateVarianceOfMean(
                                        const SampleMoment<1,T>& first_moment,
                                        const SampleMoment<2,T>& second_moment,
                                        const size_t number_of_samples );

//! Calculate the relative error
template<typename T>
double calculateRelativeError( const SampleMoment<1,T>& first_moment,
                               const SampleMoment<2,T>& second_moment,
                               const size_t number_of_samples );

//! Calculate the relative variance of the variance (VOV)
template<typename T>
double calculateRelativeVOV( const SampleMoment<1,T>& first_moment,
                             const SampleMoment<2,T>& second_moment,
                             const SampleMoment<3,T>& third_moment,
                             const SampleMoment<4,T>& fourth_moment,
                             const size_t number_of_samples );

//! Calculate the figure of merit (FOM)
double calculateFOM( const double relative_error, const double time );
  
} // end Utility namespace

//---------------------------------------------------------------------------//
// Template Includes.
//---------------------------------------------------------------------------//

#include "Utility_SampleMoment_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_SAMPLE_MOMENT_HPP

//---------------------------------------------------------------------------//
// end Utility_SampleMoment.hpp
//---------------------------------------------------------------------------//
