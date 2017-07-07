//---------------------------------------------------------------------------//
//!
//! \file   Utility_UniformDistribution_def.hpp
//! \author Alex Robinson
//! \brief  Uniform distribution class definition
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_UNIFORM_DISTRIBUTION_DEF_HPP
#define UTILITY_UNIFORM_DISTRIBUTION_DEF_HPP

// FRENSIE Includes
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_ExceptionCatchMacros.hpp"
#include "Utility_ExplicitTemplateInstantiationMacros.hpp"

namespace Utility{

// Explicit instantiation (extern declaration)
EXTERN_EXPLICIT_TEMPLATE_CLASS_INST( UnitAwareUniformDistribution<void,void> );

// Default constructor
template<typename IndependentUnit, typename DependentUnit>
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::UnitAwareUniformDistribution()
{ /* ... */ }

// Constructor
/*! \details A quantity with a different unit can be used as an input. This
 * will be explicitly cast to the desired unit during object construction.
 */
template<typename IndependentUnit, typename DependentUnit>
template<typename InputIndepQuantity, typename InputDepQuantity>
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::UnitAwareUniformDistribution(
	                       const InputIndepQuantity& min_independent_value,
			       const InputIndepQuantity& max_independent_value,
			       const InputDepQuantity& dependent_value )
  : d_min_independent_value( min_independent_value ),
    d_max_independent_value( max_independent_value ),
    d_dependent_value( dependent_value ),
    d_pdf_value()
{
  // Make sure that the values are valid
  testPrecondition( !QuantityTraits<InputIndepQuantity>::isnaninf(
						     min_independent_value ) );
  testPrecondition( !QuantityTraits<InputIndepQuantity>::isnaninf(
						     max_independent_value ) );
  testPrecondition( !QuantityTraits<InputDepQuantity>::isnaninf(
							   dependent_value ) );
  // Make sure that the max value is greater than the min value
  testPrecondition( max_independent_value > min_independent_value );

  // Calculate the pdf value
  this->calculatePDFValue();
}

// Copy constructor
/*! \details Just like boost::units::quantity objects, the unit-aware
 * distribution can be explicitly cast to a distribution with compatible
 * units. If the units are not compatible, this function will not compile. Note
 * that this allows distributions to be scaled safely (unit conversions
 * are completely taken care of by boost::units)!
 * Note: The Dummy template parameter is used to keep the Python interfaces
 *       generated by SWIG clean.
 */
template<typename IndependentUnit, typename DependentUnit>
template<typename InputIndepUnit, typename InputDepUnit, typename Dummy>
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::UnitAwareUniformDistribution( const UnitAwareUniformDistribution<InputIndepUnit,InputDepUnit>& dist_instance )
  : d_min_independent_value( dist_instance.d_min_independent_value ),
    d_max_independent_value( dist_instance.d_max_independent_value ),
    d_dependent_value( dist_instance.d_dependent_value ),
    d_pdf_value()
{
  // Make sure that the values are valid
  remember( typedef QuantityTraits<typename UnitAwareUniformDistribution<InputIndepUnit,InputDepUnit>::IndepQuantity> InputIQT );
  remember( typedef QuantityTraits<typename UnitAwareUniformDistribution<InputIndepUnit,InputDepUnit>::DepQuantity> InputDQT );
  testPrecondition( !InputIQT::isnaninf( dist_instance.d_min_independent_value ) );
  testPrecondition( !InputIQT::isnaninf( dist_instance.d_max_independent_value ) );
  testPrecondition( !InputDQT::isnaninf( dist_instance.d_dependent_value ) );
  // Make sure that the max value is greater than the min value
  testPrecondition( dist_instance.d_max_independent_value >
		    dist_instance.d_min_independent_value );

  // Calculate the pdf value
  this->calculatePDFValue();
}

// Copy constructor (copying from unitless distribution only)
template<typename IndependentUnit, typename DependentUnit>
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::UnitAwareUniformDistribution( const UnitAwareUniformDistribution<void,void>& unitless_dist_instance, int )
  : d_min_independent_value( IQT::initializeQuantity( unitless_dist_instance.d_min_independent_value ) ),
    d_max_independent_value( IQT::initializeQuantity( unitless_dist_instance.d_max_independent_value ) ),
    d_dependent_value( DQT::initializeQuantity( unitless_dist_instance.d_dependent_value ) ),
    d_pdf_value()
{
  // Make sure the values are valid
  testPrecondition( !QT::isnaninf( unitless_dist_instance.d_min_independent_value ) );
  testPrecondition( !QT::isnaninf( unitless_dist_instance.d_max_independent_value ) );
  testPrecondition( !QT::isnaninf( unitless_dist_instance.d_dependent_value ));
  // Make sure that the max value is greater than the min value
  testPrecondition( unitless_dist_instance.d_max_independent_value >
		    unitless_dist_instance.d_min_independent_value );

  // Calculate the pdf value
  this->calculatePDFValue();
}

// Construct distribution from a unitless dist. (potentially dangerous)
/*! \details Constructing a unit-aware distribution from a unitless
 * distribution is potentially dangerous. By forcing users to construct objects
 * using this method instead of a standard constructor we are trying to make
 * sure users are aware of the danger. This is designed to mimic the interface
 * of the boost::units::quantity, which also has to deal with this issue.
 */
template<typename IndependentUnit, typename DependentUnit>
UnitAwareUniformDistribution<IndependentUnit,DependentUnit> UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::fromUnitlessDistribution( const UnitAwareUniformDistribution<void,void>& unitless_distribution )
{
  return ThisType( unitless_distribution, 0 );
}

// Assignment operator
template<typename IndependentUnit, typename DependentUnit>
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>& UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::operator=(
       const UnitAwareUniformDistribution<IndependentUnit,DependentUnit>& dist_instance )
{
  // Make sure that the distribution is valid
  testPrecondition( dist_instance.d_max_independent_value >
		    dist_instance.d_min_independent_value );

  if( this != &dist_instance )
  {
    d_min_independent_value = dist_instance.d_min_independent_value;
    d_max_independent_value = dist_instance.d_max_independent_value;
    d_dependent_value = dist_instance.d_dependent_value;
    d_pdf_value = dist_instance.d_pdf_value;
  }

  return *this;
}

// Evaluate the distribution
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::DepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::evaluate(
const typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity indep_var_value ) const
{
  if( indep_var_value >= d_min_independent_value &&
      indep_var_value <= d_max_independent_value )
    return d_dependent_value;
  else
    return QuantityTraits<DepQuantity>::zero();
}

// Evaluate the PDF
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::InverseIndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::evaluatePDF(
const typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity indep_var_value ) const
{
  if( indep_var_value >= d_min_independent_value &&
      indep_var_value <= d_max_independent_value )
    return d_pdf_value;
  else
    return QuantityTraits<InverseIndepQuantity>::zero();
}

// Evaluate the CDF
template<typename IndependentUnit, typename DependentUnit>
double UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::evaluateCDF(
const typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity indep_var_value ) const
{
  if( indep_var_value >= d_min_independent_value &&
      indep_var_value <= d_max_independent_value )
    return d_pdf_value*(indep_var_value - d_min_independent_value);
  else if( indep_var_value < d_min_independent_value )
    return 0.0;
  else
    return 1.0;
}

// Return a random sample from the distribution
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::sample() const
{
  return ThisType::sample( d_min_independent_value, d_max_independent_value );
}

// Return a random sample from the distribution
template<typename IndependentUnit, typename DependentUnit>
inline typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::sample(
 const typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity min_independent_value,
 const typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity max_independent_value )
{
  // Make sure that the max value is greater than the min value
  testPrecondition( max_independent_value > min_independent_value );

  return ThisType::sampleWithRandomNumber(
			    min_independent_value,
			    max_independent_value,
			    RandomNumberGenerator::getRandomNumber<double>() );
}

// Return a random sample from the corresponding CDF and record the number of trials
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::sampleAndRecordTrials( DistributionTraits::Counter& trials ) const
{
  return ThisType::sampleAndRecordTrials( d_min_independent_value,
					  d_max_independent_value,
					  trials );
}

// Return a random sample from the distribution and record the number of trials
template<typename IndependentUnit, typename DependentUnit>
inline typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::sampleAndRecordTrials(
 const typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity min_independent_value,
 const typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity max_independent_value,
 DistributionTraits::Counter& trials )
{
  // Make sure that the max value is greater than the min value
  testPrecondition( max_independent_value > min_independent_value );

  ++trials;

  return ThisType::sample( min_independent_value, max_independent_value );
}

// Return a random sample from the distribution at the given CDF value
template<typename IndependentUnit, typename DependentUnit>
inline typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::sampleWithRandomNumber(
					     const double random_number ) const
{
  // Make sure the random number is valid
  testPrecondition( random_number >= 0.0 );
  testPrecondition( random_number <= 1.0 );

  return ThisType::sampleWithRandomNumber( d_min_independent_value,
					   d_max_independent_value,
					   random_number );
}

// Return a random sample from the distribution at the given CDF value
template<typename IndependentUnit, typename DependentUnit>
inline typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::sampleWithRandomNumber(
  const typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity min_independent_value,
  const typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity max_independent_value,
  const double random_number )
{
  // Make sure that the max value is greater than the min value
  testPrecondition( max_independent_value > min_independent_value );
  // Make sure the random number is valid
  testPrecondition( random_number >= 0.0 );
  testPrecondition( random_number <= 1.0 );

  return random_number*(max_independent_value - min_independent_value) +
    min_independent_value;
}

// Return a random sample from the distribution at the given CDF value in a subrange
template<typename IndependentUnit, typename DependentUnit>
inline typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::sampleWithRandomNumberInSubrange(
	    const double random_number,
	    const typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity max_indep_var ) const
{
  // Make sure the random number is valid
  testPrecondition( random_number >= 0.0 );
  testPrecondition( random_number <= 1.0 );
  // Make sure the upper bound of the subrange is valid
  testPrecondition( max_indep_var <= d_max_independent_value );
  testPrecondition( max_indep_var >= d_min_independent_value );

  return ThisType::sampleWithRandomNumber( d_min_independent_value,
					   max_indep_var,
					   random_number );
}

// Return a random sample and sampled index from the corresponding CDF
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::sampleAndRecordBinIndex(
					    unsigned& sampled_bin_index ) const
{
  sampled_bin_index = 0u;

  return ThisType::sample( d_min_independent_value, d_max_independent_value );
}

// Return a random sample from the corresponding CDF in a subrange
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::sampleInSubrange(const UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity max_indep_var ) const
{
  // Make sure the upper bound of the subrange is valid
  testPrecondition( max_indep_var <= d_max_independent_value );
  testPrecondition( max_indep_var >= d_min_independent_value );

  return ThisType::sampleWithRandomNumber(
			    d_min_independent_value,
			    max_indep_var,
			    RandomNumberGenerator::getRandomNumber<double>() );
}

// Return the upper bound of the distribution independent variable
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::getUpperBoundOfIndepVar() const
{
  return d_max_independent_value;
}

// Return the lower bound of the distribution independent variable
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::getLowerBoundOfIndepVar() const
{
  return d_min_independent_value;
}

// Return the distribution type
template<typename IndependentUnit, typename DependentUnit>
OneDDistributionType UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::getDistributionType() const
{
  return UnitAwareUniformDistribution::distribution_type;
}

// Test if the distribution is continuous
template<typename IndependentUnit, typename DependentUnit>
bool UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::isContinuous() const
{
  return true;
}

// Method for placing the object in an output stream
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::toStream( std::ostream& os ) const
{
  os << "{" << Utility::convertOneDDistributionTypeToString( UnitAwareUniformDistribution::distribution_type ) << ", ";

  Utility::toStream( os, getRawQuantity( d_min_independent_value ) );

  os << ", ";

  Utility::toStream( os, getRawQuantity( d_max_independent_value ) );
  
  os << ", ";

  Utility::toStream( os, getRawQuantity( d_dependent_value ) );
  
  os << "}";
}

// Method for initializing the object from an input stream
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::fromStream(
                                                           std::istream& is,
                                                           const std::string& )
{
  std::string dist_rep;
  std::getline( is, dist_rep, '}' );
  dist_rep += '}';

  VariantVector distribution_data;

  try{
    distribution_data =
      Utility::variant_cast<VariantVector>( Utility::Variant::fromValue( dist_rep ) );
  }
  EXCEPTION_CATCH_RETHROW( Utility::StringConversionException,
                           "Could not extract the distribution data from "
                           "the stream!" );

  // Verify that the correct amount of distribution data is present
  TEST_FOR_EXCEPTION( distribution_data.size() < 3 ||
                      distribution_data.size() > 4,
                      Utility::StringConversionException,
                      "The uniform distribution cannot be constructed "
                      "because the string representation is not valid!" );

  // Verify that the distribution type is uniform
  this->verifyDistributionType( distribution_data[0] );

  // Extract the min independent value
  this->setMinIndependentValue( distribution_data[1] );

  // Extract the max independent value
  this->setMaxIndependentValue( distribution_data[2] );

  // Verify that the independent values are valid
  this->verifyValidIndependentValues();

  // Extract the dependent value
  if( distribution_data.size() == 4 )
    this->setDependentValue( distribution_data[3] );
  else
    this->setDependentValue( Utility::Variant::fromValue( 1.0 ) );

  // Calculate the distribution's pdf
  this->calculatePDFValue();
}

// Method for placing the object in the desired property tree node
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::toNode(
                                                 const std::string& node_key,
                                                 Utility::PropertyTree& ptree,
                                                 const bool inline_data ) const
{
  if( inline_data )
  {
    std::ostringstream oss;

    this->toStream( oss );

    ptree.put( node_key, oss.str() );
  }
  else
  {
    Utility::PropertyTree child_tree;
  
    child_tree.put( "type", Utility::convertOneDDistributionTypeToString( UnitAwareUniformDistribution::distribution_type ) );
    child_tree.put( "min indep value", getRawQuantity( d_min_independent_value ) );
    child_tree.put( "max indep value", getRawQuantity( d_max_independent_value ) );
    child_tree.put( "dep value", getRawQuantity( d_dependent_value ) );
    
    ptree.put_child( node_key, child_tree );
  }
}

// Method for initializing the object from a property tree node
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::fromNode(
                                    const Utility::PropertyTree& node,
                                    std::vector<std::string>& unused_children )
{
  // Initialize from inline data
  if( node.size() == 0 )
  {
    std::istringstream iss( node.data().toString() );

    this->fromStream( iss );
  }
  // Initialize from child nodes
  else
  {
    Utility::PropertyTree::const_iterator node_it, node_end;
    node_it = node.begin();
    node_end = node.end();
    
    bool type_verified = false;
    bool min_indep_val_set = false;
    bool max_indep_val_set = false;
    bool dep_val_set = false;
    
    while( node_it != node_end )
    {
      std::string child_node_key =
        boost::algorithm::to_lower_copy( node_it->first );

      // Verify the distribution type
      if( child_node_key.find( "type" ) < child_node_key.size() )
      {
        this->verifyDistributionType( node_it->second.data() );
        
        type_verified = true;
      }
      
      // Extract the min indep value
      else if( child_node_key.find( "min" ) < child_node_key.size() )
      {
        this->setMinIndependentValue( node_it->second.data() );
        
        min_indep_val_set = true;
      }

      // Extract the max indep value
      else if( child_node_key.find( "max" ) < child_node_key.size() )
      {
        this->setMaxIndependentValue( node_it->second.data() );
        
        max_indep_val_set = true;
      }
      
      // Extract the dependent value
      else if( child_node_key.find( "dep" ) < child_node_key.size() )
      {
        this->setDependentValue( node_it->second.data() );
        
        dep_val_set = true;
      }

      // This child node is unused (and is not a comment)
      else if( child_node_key.find( PTREE_COMMENT_NODE_KEY ) >=
               child_node_key.size() )
      {
        unused_children.push_back( node_it->first );
      }
      
      ++node_it;
    }

    // Make sure that the distribution type was verified
    TEST_FOR_EXCEPTION( !type_verified,
                        std::runtime_error,
                        "The uniform distribution could not be constructed "
                        "because the type could not be verified!" );
    
    // Make sure that the min indep value was set
    TEST_FOR_EXCEPTION( !min_indep_val_set,
                        std::runtime_error,
                        "The uniform distribution could not be constructed "
                        "because the min independent value was not "
                        "specified!" );
    
    // Make sure that the max indep value was set
    TEST_FOR_EXCEPTION( !max_indep_val_set,
                        std::runtime_error,
                        "The uniform distribution could not be constructed "
                        "because the max independent value was not "
                        "specified!" );
    
    // Verify that the independent values are valid
    this->verifyValidIndependentValues();
    
    // Check if the dependent value was set
    if( !dep_val_set )
      this->setDependentValue( Utility::Variant::fromValue( 1.0 ) );

    // Calculate the distribution's pdf
    this->calculatePDFValue();
  }
}

// Verify the distribution type
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::verifyDistributionType( const Utility::Variant& type_data )
{
  std::string distribution_type = type_data.toString();
  boost::algorithm::to_lower( distribution_type );
  
  TEST_FOR_EXCEPTION( distribution_type.find( "uniform" ) >=
                      distribution_type.size(),
                      Utility::StringConversionException,
                      "The uniform distribution cannot be constructed "
                      "because the distribution type ("
                      << distribution_type << ") does not match!" );
}

// Set the min indep value
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::setMinIndependentValue( const Utility::Variant& min_indep_data )
{
  try{
    Utility::setQuantity( d_min_independent_value,
                          Utility::variant_cast<double>( min_indep_data ) );
  }
  EXCEPTION_CATCH_RETHROW( Utility::StringConversionException,
                           "Could not extract the min independent value!" );

  // Verify that the min independent value is valid
  TEST_FOR_EXCEPTION( IQT::isnaninf( d_min_independent_value ),
		      InvalidDistributionStringRepresentation,
		      "The uniform distribution cannot be "
		      "constructed because of an invalid min "
		      "independent value " <<
		      d_min_independent_value );
}

// Set the max indep value
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::setMaxIndependentValue( const Utility::Variant& max_indep_data )
{
  try{
    Utility::setQuantity( d_max_independent_value,
                          Utility::variant_cast<double>( max_indep_data ) );
  }
  EXCEPTION_CATCH_RETHROW( Utility::StringConversionException,
                           "Could not extract the max independent value!" );

  // Verify that the max independent value is valid
  TEST_FOR_EXCEPTION( IQT::isnaninf( d_max_independent_value ),
		      Utility::StringConversionException,
		      "The uniform distribution cannot be "
		      "constructed because of an invalid max "
		      "independent value " <<
		      d_max_independent_value );
}

// Verify that the independent values are valid
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::verifyValidIndependentValues()
{
  TEST_FOR_EXCEPTION( d_max_independent_value <= d_min_independent_value,
		      Utility::StringConversionException,
		      "The uniform distribution cannot be constructed because "
                      "of invalid independent values!" );
}

// Set the dependent indep value
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::setDependentValue( const Utility::Variant& dep_data )
{
  try{
    Utility::setQuantity( d_dependent_value,
                          Utility::variant_cast<double>( dep_data ) );
  }
  EXCEPTION_CATCH_RETHROW( Utility::StringConversionException,
                           "Could not extract the dependent value!" );

  TEST_FOR_EXCEPTION( DQT::isnaninf( d_dependent_value ),
		      Utility::StringConversionException,
		      "The uniform distribution cannot be "
		      "constructed because of an invalid dependent "
		      "value " << d_dependent_value );
}

// Calculate the PDF value
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::calculatePDFValue()
{
  d_pdf_value = 1.0/(d_max_independent_value - d_min_independent_value);
}

// Test if the dependent variable can be zero within the indep bounds
template<typename IndependentUnit, typename DependentUnit>
bool UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::canDepVarBeZeroInIndepBounds() const
{
  return false;
}

// Equality comparison operator
template<typename IndependentUnit, typename DependentUnit>
bool UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::operator==( const UnitAwareUniformDistribution& other ) const
{
  return d_min_independent_value == other.d_min_independent_value &&
    d_max_independent_value == other.d_max_independent_value &&
    d_dependent_value == other.d_dependent_value;
}

// Inequality comparison operator
template<typename IndependentUnit, typename DependentUnit>
bool UnitAwareUniformDistribution<IndependentUnit,DependentUnit>::operator!=( const UnitAwareUniformDistribution& other ) const
{
  return d_min_independent_value != other.d_min_independent_value ||
    d_max_independent_value != other.d_max_independent_value ||
    d_dependent_value != other.d_dependent_value;
}

} // end Utility namespace

#endif // end UTILITY_UNIFORM_DISTRIBUTION_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_UniformDistribution_def.hpp
//---------------------------------------------------------------------------//
