//---------------------------------------------------------------------------//
//!
//! \file   Utility_DeltaDistribution_def.hpp
//! \author Alex Robinson
//! \brief  Delta distribution class declaration.
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_DELTA_DISTRIBUTION_DEF_HPP
#define UTILITY_DELTA_DISTRIBUTION_DEF_HPP

// Std Lib Includes
#include <limits>

// FRENSIE Includes
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_ExceptionCatchMacros.hpp"
#include "Utility_ExplicitTemplateInstantiationMacros.hpp"
#include "Utility_ContractException.hpp"

namespace Utility{

// Explicit instantiation (extern declaration)
EXTERN_EXPLICIT_TEMPLATE_CLASS_INST( UnitAwareDeltaDistribution<void,void> );

// Default constructor
template<typename IndependentUnit, typename DependentUnit>
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::UnitAwareDeltaDistribution()
{ /* ... */ }

// Basic Constructor
template<typename IndependentUnit, typename DependentUnit>
template<typename InputIndepQuantity>
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::UnitAwareDeltaDistribution( const InputIndepQuantity location )
  : d_location( location ),
    d_multiplier( DQT::one() )
{
  // Make sure that the point is valid
  testPrecondition( !QuantityTraits<InputIndepQuantity>::isnaninf(location) );
}

// Advanced Constructor
template<typename IndependentUnit, typename DependentUnit>
template<typename InputIndepQuantity, typename InputDepQuantity>
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::UnitAwareDeltaDistribution(
                                            const InputIndepQuantity location,
					    const InputDepQuantity multiplier )
  : d_location( location ),
    d_multiplier( multiplier )
{
  // Make sure that the point is valid
  testPrecondition( !QuantityTraits<InputIndepQuantity>::isnaninf( location ));
  // Make sure the multiplier is valid
  testPrecondition( !QuantityTraits<InputDepQuantity>::isnaninf( multiplier ));
  testPrecondition( multiplier != QuantityTraits<InputDepQuantity>::zero() );
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
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::UnitAwareDeltaDistribution(
 const UnitAwareDeltaDistribution<InputIndepUnit,InputDepUnit>& dist_instance )
  : d_location( dist_instance.d_location ),
    d_multiplier( dist_instance.d_multiplier )
{
  remember( typedef QuantityTraits<typename UnitAwareDeltaDistribution<InputIndepUnit,InputDepUnit>::IndepQuantity> InputIQT );
  remember( typedef QuantityTraits<typename UnitAwareDeltaDistribution<InputIndepUnit,InputDepUnit>::DepQuantity> InputDQT );
  // Make sure that the point is valid
  testPrecondition( !InputIQT::isnaninf( dist_instance.d_location ) );
  // Make sure that the multiplier is valid
  testPrecondition( !InputDQT::isnaninf( dist_instance.d_multiplier ) );
  testPrecondition( dist_instance.d_multiplier != InputDQT::zero() );
}

// Copy constructor (copying from unitless distribution only)
template<typename IndependentUnit, typename DependentUnit>
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::UnitAwareDeltaDistribution(
const UnitAwareDeltaDistribution<void,void>& unitless_dist_instance, int )
  : d_location( IQT::initializeQuantity( unitless_dist_instance.d_location ) ),
    d_multiplier( DQT::initializeQuantity( unitless_dist_instance.d_multiplier ) )
{
  // Make sure that the point is valid
  testPrecondition( !QT::isnaninf( unitless_dist_instance.d_location ) );
  // Make sure that the multiplier is valid
  testPrecondition( !QT::isnaninf( unitless_dist_instance.d_multiplier ) );
  testPrecondition( unitless_dist_instance.d_multiplier != 0.0 );
}

// Construct distribution from a unitless dist. (potentially dangerous)
/*! \details Constructing a unit-aware distribution from a unitless
 * distribution is potentially dangerous. By forcing users to construct objects
 * using this method instead of a standard constructor we are trying to make
 * sure users are aware of the danger. This is designed to mimic the interface
 * of the boost::units::quantity, which also has to deal with this issue.
 */
template<typename IndependentUnit, typename DependentUnit>
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit> UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::fromUnitlessDistribution(
	   const UnitAwareDeltaDistribution<void,void>& unitless_distribution )
{
  return ThisType( unitless_distribution, 0 );
}

// Assignment operator
template<typename IndependentUnit, typename DependentUnit>
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>&
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::operator=(
 const UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>& dist_instance )
{
  if( this != &dist_instance )
  {
    d_location = dist_instance.d_location;

    d_multiplier = dist_instance.d_multiplier;
  }

  return *this;
}

// Evaluate the distribution
/*! \details The delta distribution is a continuous distribution
 * that can evaluate to one of two values: 0.0 and infinity. It is more
 * useful (and safer) to simply return 0.0 and a (the multiplier).
 */
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::DepQuantity
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::evaluate(
   const typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity indep_var_value ) const
{
  if( indep_var_value == d_location )
    return d_multiplier;
  else
    return DQT::zero();
}

// Evaluate the PDF
/*! \details The delta distribution is a continuous distribution
 * that can evaluate to one of two values: 0.0 and infinity. It is more
 * useful (and safer) to simply return 0.0 and 1.0.
 */
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::InverseIndepQuantity
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::evaluatePDF(
 const typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity indep_var_value ) const
{
  if( indep_var_value == d_location )
    return QuantityTraits<InverseIndepQuantity>::one();
  else
    return QuantityTraits<InverseIndepQuantity>::zero();
}

// Evaluate the CDF
template<typename IndependentUnit, typename DependentUnit>
double UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::evaluateCDF(
  const typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity indep_var_value ) const
{
  if( indep_var_value < d_location )
    return 0.0;
  else
    return 1.0;
}

// Return a random sample from the distribution
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::sample() const
{
  return d_location;
}

// Return a random sample from the corresponding CDF and record the number of trials
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::sampleAndRecordTrials( DistributionTraits::Counter& trials ) const
{
  ++trials;

  return d_location;
}

// Return a random sample from the distribution and the sampled index
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::sampleAndRecordBinIndex(
					    unsigned& sampled_bin_index ) const
{
  sampled_bin_index = 0;

  return d_location;
}

// Return a random sample from the distribution at the given CDF value
/*! \details The random number will be ignored since only a single value can
 * every be returned
 */
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::sampleWithRandomNumber(
					     const double random_number ) const
{
  return d_location;
}

// Return a random sample from the distribution in a subrange
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::sampleInSubrange(
 const typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity max_indep_var ) const
{
  // Make sure the max independent variable is valid
  testPrecondition( max_indep_var >= d_location );

  return d_location;
}

// Return a random sample from the distribution at the given CDF value in a subrange
/*! \details The random number will be ignored since only a single value can
 * ever be returned
 */
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::sampleWithRandomNumberInSubrange(
 const double random_number,
 const typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity max_indep_var ) const
{
  // Make sure the max independent variable is valid
  testPrecondition( max_indep_var >= d_location );

  return d_location;
}

// Return the maximum point at which the distribution is non-zero
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::getUpperBoundOfIndepVar() const
{
  return d_location;
}

// Return the minimum point at which the distribution is non-zero
template<typename IndependentUnit, typename DependentUnit>
typename UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::IndepQuantity
UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::getLowerBoundOfIndepVar() const
{
  return d_location;
}

// Return the distribution type
template<typename IndependentUnit, typename DependentUnit>
OneDDistributionType UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::getDistributionType() const
{
  return ThisType::distribution_type;
}

// Return the distribution type name
template<typename IndependentUnit, typename DependentUnit>
std::string UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::getDistributionTypeName(
                                                       const bool verbose_name,
                                                       const bool lowercase )
{
  std::string name = "Delta";

  if( verbose_name )
    name += " Distribution";

  if( lowercase )
    boost::algorithm::to_lower( name );

  return name;
}

// Check if the type name matches the distribution type name
/*! \detail The type name comparison is case-insensitive. A positive match
 * will be reported if the type name has a substring equal to "delta".
 */
template<typename IndependentUnit, typename DependentUnit>
bool UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::doesTypeNameMatch( const std::string type_name )
{
  std::string lower_type_name = boost::algorithm::to_lower_copy( type_name );
  
  return lower_type_name.find(ThisType::getDistributionTypeName( false, true )) < lower_type_name.size();
}

// Test if the distribution is continuous
/*! \details Though the delta distribution is technically continuous
 * because it is only non-zero at the specified point it will be treated as
 * a discrete distribution.
 */
template<typename IndependentUnit, typename DependentUnit>
bool UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::isContinuous() const
{
  return false;
}

// Method for placing the object in an output stream
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::toStream( std::ostream& os ) const
{
  os << Utility::container_start_char
     << Utility::convertOneDDistributionTypeToString( ThisType::distribution_type )
     << Utility::next_container_element_char << " ";

  Utility::toStream( os, Utility::getRawQuantity( d_location ) );

  if( d_multiplier != DQT::one() )
  {
    os << Utility::next_container_element_char << " ";
    
    Utility::toStream( os, Utility::getRawQuantity( d_multiplier ) );
  }

  os << Utility::container_end_char;
}

// Method for initializing the object from an input stream
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::fromStream(
                                                           std::istream& is,
                                                           const std::string& )
{
  // Read in the distribution representation
  std::string dist_rep;
  std::getline( is, dist_rep, Utility::container_end_char );
  dist_rep += Utility::container_end_char;

  VariantVector distribution_data;

  try{
    distribution_data =
      Utility::variant_cast<VariantVector>( Utility::Variant::fromValue( dist_rep ) );
  }
  EXCEPTION_CATCH_RETHROW( Utility::StringConversionException,
                           "Could not extract the distribution data from "
                           "the stream!" );

  // Verify that the correct amount of distribution data is present
  TEST_FOR_EXCEPTION( distribution_data.size() < 2 ||
                      distribution_data.size() > 3,
		      Utility::StringConversionException,
		      "The delta distribution cannot be constructed "
		      "because the string representation is not valid!" );

  // Verify that the distribution type is delta
  this->verifyDistributionType( distribution_data[0] );

  // Extract the location value
  this->setLocationValue( distribution_data[1] );

  // Extract the multiplier value
  if( distribution_data.size() > 2 )
    this->setMultiplierValue( distribution_data[2] );
  else
    d_multiplier = DQT::one();
}

// Method for placing the object in the desired property tree node
template<typename IndependentUnit, typename DependentUnit>
Utility::PropertyTree UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::toPropertyTree(
                                                 const bool inline_data ) const
{
  Utility::PropertyTree ptree;
  
  if( inline_data )
    ptree.put_value( *this );
  else
  {
    ptree.put( "type", Utility::convertOneDDistributionTypeToString( ThisType::distribution_type ) );
    ptree.put( "location", Utility::getRawQuantity( d_location ) );

    if( d_multiplier != DQT::one() )
      ptree.put( "multiplier", Utility::getRawQuantity( d_multiplier ) );
  }

  return ptree;
}

// Method for initializing the object from a property tree
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::fromPropertyTree(
                                    const Utility::PropertyTree& node,
                                    std::vector<std::string>& unused_children )
{
  // Initialize from inline data
  if( node.size() == 0 )
  {
    std::istringstream iss( node.data().toString() );

    try{
      this->fromStream( iss );
    }
    EXCEPTION_CATCH_RETHROW_AS( std::runtime_error,
                                Utility::PropertyTreeConversionException,
                                "Could not create the delta "
                                "distribution!" );
  }
  // Initialize from child nodes
  else
  {
    Utility::PropertyTree::const_iterator node_it, node_end;
    node_it = node.begin();
    node_end = node.end();

    bool type_verified = false;
    bool location_set = false;
    bool multiplier_set = false;

    while( node_it != node_end )
    {
      std::string child_node_key =
        boost::algorithm::to_lower_copy( node_it->first );

      // Verify the distribution type
      if( child_node_key.find( "type" ) < child_node_key.size() )
      {
        try{
          this->verifyDistributionType( node_it->second.data() );
        }
        EXCEPTION_CATCH_RETHROW_AS( std::runtime_error,
                                    Utility::PropertyTreeConversionException,
                                    "Could not create the delta "
                                    "distribution!" );

        type_verified = true;
      }

      // Extract the location value
      else if( child_node_key.find( "loc" ) < child_node_key.size() )
      {
        try{
          this->setLocationValue( node_it->second.data() );
        }
        EXCEPTION_CATCH_RETHROW_AS( std::runtime_error,
                                    Utility::PropertyTreeConversionException,
                                    "Could not create the delta "
                                    "distribution!" );

        location_set = true;
      }

      // Extract the multiplier value
      else if( child_node_key.find( "mult" ) < child_node_key.size() )
      {
        try{
          this->setMultiplierValue( node_it->second.data() );
        }
        EXCEPTION_CATCH_RETHROW_AS( std::runtime_error,
                                    Utility::PropertyTreeConversionException,
                                    "Could not create the delta "
                                    "distribution!" );

        multiplier_set = true;
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
                        Utility::PropertyTreeConversionException,
                        "The delta distribution could not be constructed "
                        "because the type could not be verified!" );

    // Make sure that the location value was set
    TEST_FOR_EXCEPTION( !location_set,
                        Utility::PropertyTreeConversionException,
                        "The delta distribution could not be constructed "
                        "because the location value was not specified!" );
    
    // Check if the multiplier was set
    if( !multiplier_set )
      d_multiplier = DQT::one();
  }
}

// Verify the distribution type
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::verifyDistributionType( const Utility::Variant& type_data )
{
  std::string distribution_type = type_data.toLowerString();

  TEST_FOR_EXCEPTION( !ThisType::doesTypeNameMatch( distribution_type ),
                      Utility::StringConversionException,
                      "The delta distribution cannot be constructed "
                      "because the distribution type ("
                      << distribution_type << ") does not match!" );
}

// Set the location value
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::setLocationValue( const Utility::Variant& location_data )
{
  try{
    Utility::setQuantity( d_location,
                          Utility::variant_cast<double>( location_data ) );
  }
  EXCEPTION_CATCH_RETHROW( Utility::StringConversionException,
                           "Could not extract the location value!" );

  // Verify that the location value is valid
  TEST_FOR_EXCEPTION( IQT::isnaninf( d_location ),
		      Utility::StringConversionException,
		      "The delta distribution cannot be constructed "
		      "because of an invalid location (" << d_location <<
		      ")!" );
}

// Set the multiplier value
template<typename IndependentUnit, typename DependentUnit>
void UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::setMultiplierValue( const Utility::Variant& multiplier_data )
{
  try{
    Utility::setQuantity( d_multiplier,
                          Utility::variant_cast<double>( multiplier_data ) );
  }
  EXCEPTION_CATCH_RETHROW( Utility::StringConversionException,
                           "Could not extract the multiplier value!" );

  TEST_FOR_EXCEPTION( DQT::isnaninf( d_multiplier ),
                      Utility::StringConversionException,
                      "The delta distribution cannot be constructed "
                      "because of an invalid multiplier ("
                      << d_multiplier << ")!" );

  TEST_FOR_EXCEPTION( d_multiplier == DQT::zero(),
                      Utility::StringConversionException,
                      "The delta distribution cannot be constructed "
                      "because of an invalid multiplier ("
                      << d_multiplier << ")!" );
}

// Equality comparison operator
template<typename IndependentUnit, typename DependentUnit>
bool UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::operator==(
 const UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>& other ) const
{
  return d_location == other.d_location &&
    d_multiplier == other.d_multiplier;
}

// Inequality comparison operator
template<typename IndependentUnit, typename DependentUnit>
bool UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::operator!=(
 const UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>& other ) const
{
  return d_location != other.d_location ||
    d_multiplier != other.d_multiplier;
}

// Test if the dependent variable can be zero within the indep bounds
template<typename IndependentUnit, typename DependentUnit>
bool UnitAwareDeltaDistribution<IndependentUnit,DependentUnit>::canDepVarBeZeroInIndepBounds() const
{
  return false;
}

} // end Utility namespace

#endif // end UTILITY_DELTA_DISTRIBUTION_DEF_HPP

//---------------------------------------------------------------------------//
// end Utiliyt_DeltaDistribution_def.hpp
//---------------------------------------------------------------------------//
