//---------------------------------------------------------------------------//
//!
//! \file   Utility_GaussKronrodIntegrator.hpp
//! \author Luke Kersting
//! \brief  Gauss-Kronrod integrator
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_GAUSS_KRONROD_INTEGRATOR_HPP
#define UTILITY_GAUSS_KRONROD_INTEGRATOR_HPP

// Std Includes
#include <queue>

// Trilinos Includes
#include <Teuchos_Array.hpp>

// FRENSIE Includes
#include "Utility_QuadratureBin.hpp"
#include "Utility_ExtrapolatedQuadratureBin.hpp"

namespace Utility{

//! The Gauss-Kronrod integrator
template<typename T>
class GaussKronrodIntegrator
{

private:

  // The quantity traits for type T
  typedef Utility::QuantityTraits<T> QT;

  // The bin traits priority queue
  typedef std::priority_queue<QuadratureBin<T> > BinQueue;

  // The bin traits array
  typedef Teuchos::Array<ExtrapolatedQuadratureBin<T> > BinArray;

public:

  //! Constructor
  GaussKronrodIntegrator( const T relative_error_tol,
                          const T absolute_error_tol = QT::zero(),
                          const size_t subinterval_limit = 1000 );

  //! Destructor
  ~GaussKronrodIntegrator()
  { /* ... */ }

  //! Integrate the function
  // template<typename Functor>
  // void integrate( Functor& integrand, 
  //       	  T lower_limit, 
  //       	  T upper_limit,
  //       	  T& result,
  //       	  T& absolute_error,
  //       	  size_t& number_of_function_evals ) const;

  //! Integrate the function adaptively with BinQueue
  template<int Points, typename Functor>
  void integrateAdaptively( Functor& integrand,
			    T lower_limit,
			    T upper_limit,
			    T& result,
			    T& absolute_error ) const;

  //! Integrate the function with point rule
  template<int Points, typename Functor>
  void integrateWithPointRule( 
                Functor& integrand,
                T lower_limit,
                T upper_limit,
                T& result,
                T& absolute_error,
                T& result_abs, 
                T& result_asc ) const;

  //! Integrate the function over a semi-infinite interval (+infinity)
  // template<typename Functor>
  // void integrateSemiInfiniteIntervalUpper( Functor& integrand,
  //       				   T lower_limit,
  //       				   T& result,
  //       				   T& absolute_error ) const;

  //! Integrate the function over a semi-infinite interval (-infinity)
  // template<typename Functor>
  // void integrateSemiInfiniteIntervalLower( Functor& integrand,
  //       				   T upper_limit,
  //       				   T& result,
  //       				   T& absolute_error ) const;

  //! Integrate the function over an infinite interval (-infinity,+infinity)
  // template<typename Functor>
  // void integrateInfiniteInterval( Functor& integrand,
  //       			  T& result,
  //       			  T& absolute_error ) const;

  //! Integrate a function with integrable singularities adaptively
  // template<typename Functor>
  // void integrateAdaptivelyWynnEpsilon( Functor& integrand,
  //       			       T lower_limit,
  //       			       T upper_limit,
  //       			       T& result,
  //       			       T& absolute_error ) const;

  //! Integrate a function with known integrable singularities adaptively
  template<typename Functor>
  void integrateAdaptivelyWynnEpsilon( 
			  Functor& integrand,
			  const Teuchos::ArrayView<T>& points_of_interest,
			  T& result,
			  T& absolute_error ) const;

protected:

  // Conduct the first iteration of the adapative integration
  template<int Points, typename Functor>
  bool integrateAdaptivelyInitialIteration( Functor& integrand,
                                            QuadratureBin<T>& bin ) const;

  // Conduct the other iterations of the adaptive integration
  template<int Points, typename Functor>
  void integrateAdaptivelyIterate( Functor& integrand,
                                   QuadratureBin<T>& initial_bin,
                                   T& result,
                                   T& absolute_error ) const;

  // Initialize the bins for the adaptive Wynn Epsilon integration
  template<int Points, typename Functor>
  bool initializeBinsWynnEpsilon(
                               Functor& integrand,
                               const Teuchos::ArrayView<T>& points_of_interest,
                               Teuchos::Array<unsigned>& bin_order,
                               T& result,
                               T& absolute_error ) const;
  
  // Calculate the quadrature upper and lower integrand values at an abscissa
  template<typename Functor>
  void calculateQuadratureIntegrandValuesAtAbscissa( 
    Functor& integrand, 
    T abscissa,
    T half_length,
    T midpoint,
    T& integrand_value_lower,
    T& integrand_value_upper ) const;

  // Bisect and integrate the given bin interval
  template<int Points, typename Functor, typename Bin>
  void bisectAndIntegrateBinInterval( 
    Functor& integrand, 
    const Bin& bin,
    Bin& bin_1,
    Bin& bin_2,
    T& bin_1_asc,
    T& bin_2_asc ) const;

  // Rescale absolute error from integration
  void rescaleAbsoluteError( 
    T& absolute_error, 
    T result_abs, 
    T result_asc ) const;

  // Test if subinterval is too small
  template<int Points>
  bool subintervalTooSmall( T& lower_limit_1, 
                            T& lower_limit_2, 
                            T& upper_limit_2 ) const;

  // check the roundoff error 
  void checkRoundoffError( 
                       const QuadratureBin<T>& bin, 
                       const QuadratureBin<T>& bin_1, 
                       const QuadratureBin<T>& bin_2,    
                       const T& bin_1_asc,
                       const T& bin_2_asc,
                       int& round_off_1,
                       int& round_off_2,
                       const int number_of_iterations ) const;

  // check the roundoff error 
  void checkRoundoffError( 
                       const ExtrapolatedQuadratureBin<T>& bin, 
                       const ExtrapolatedQuadratureBin<T>& bin_1, 
                       const ExtrapolatedQuadratureBin<T>& bin_2,    
                       const T& bin_1_asc,
                       const T& bin_2_asc,
                       int& round_off_1,
                       int& round_off_2,
                       int& round_off_3,
                       const bool extrapolate, 
                       const int number_of_iterations ) const;
 
  // Sort the bin order from highest to lowest error 
  void sortBins( 
        Teuchos::Array<int>& bin_order,
        BinArray& bin_array, 
        const ExtrapolatedQuadratureBin<T>& bin_1,
        const ExtrapolatedQuadratureBin<T>& bin_2,
        const int& number_of_intervals,
        int& nr_max ) const;

  // get the Wynn Epsilon-Algoirithm extrapolated value
  void getWynnEpsilonAlgorithmExtrapolation( 
        Teuchos::Array<T>& bin_extrapolated_result, 
        Teuchos::Array<T>& last_three_results, 
        T& extrapolated_result, 
        T& extrapolated_error,  
        int& number_of_extrapolated_intervals,
        int& number_of_extrapolated_calls  ) const;

private:
  // The relative error tolerance
  T d_relative_error_tol;

  // The absolute error tolerance
  T d_absolute_error_tol;

  // The subinterval limit
  size_t d_subinterval_limit;

  // return epsilon numerical limit for type T
  T getLimitEpsilon() const;

  // return min numerical limit for type T
  T getLimitMin() const;

  // return max numerical limit for type T
  T getLimitMax() const;

  // return max of two variables of type T
  T getMax( T variable_1, T variable_2 ) const;
};

} // end Utility namespace

//---------------------------------------------------------------------------//
// Template Includes
//---------------------------------------------------------------------------//

#include "Utility_GaussKronrodIntegrator_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_GUASS_KRONROD_INTEGRATOR_HPP

//---------------------------------------------------------------------------//
// end Utility_GaussKronrodIntegrator.hpp
//---------------------------------------------------------------------------//
