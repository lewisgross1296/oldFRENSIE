//---------------------------------------------------------------------------//
//!
//! \file   Utility_BivariateDistributionHelpers.i
//! \author Luke Kersting
//! \brief  The distribution helper macros
//!
//---------------------------------------------------------------------------//

// Include std::string support
%include <std_string.i>

 // Import the PyFrensie_ArraySharedPtr.i
%include "PyFrensie_ArraySharedPtr.i"

//---------------------------------------------------------------------------//
// Helper macros for getting the distribution name
//---------------------------------------------------------------------------//
// Get the name of a distribution that has 3,4,5,6 template parameters
// It must be included in both the interface file and the wrapper file since
// the preprocessor will have to parse it separately in both files
%{
#define _BI_DIST_NAME_6_ARGS( dist_base_name, arg_1, arg_2, arg_3, arg_4, arg_5, arg_6 ) Utility::UnitAware ## dist_base_name<arg_2<arg_1>, arg_3<arg_4,arg_5,arg_6> >
%}
#define _BI_DIST_NAME_6_ARGS( dist_base_name, arg_1, arg_2, arg_3, arg_4, arg_5, arg_6 ) Utility::UnitAware ## dist_base_name<arg_2<arg_1>, arg_3<arg_4,arg_5,arg_6> >

%{
#define _BI_DIST_NAME_5_ARGS( dist_base_name, arg_1, arg_2, arg_3, arg_4, arg_5 ) Utility::UnitAware ## dist_base_name<arg_2<arg_1>,arg_3,arg_4,arg_5>
%}
#define _BI_DIST_NAME_5_ARGS( dist_base_name, arg_1, arg_2, arg_3, arg_4, arg_5 ) Utility::UnitAware ## dist_base_name<arg_2<arg_1>,arg_3,arg_4,arg_5>

%{
#define _BI_DIST_NAME_4_ARGS( dist_base_name, arg_1, arg_2, arg_3, arg_4 ) Utility::UnitAware ## dist_base_name<arg_1,arg_2,arg_3,arg_4>
%}
#define _BI_DIST_NAME_4_ARGS( dist_base_name, arg_1, arg_2, arg_3, arg_4 ) Utility::UnitAware ## dist_base_name<arg_1,arg_2,arg_3,arg_4>

%{
#define _BI_DIST_NAME_3_ARGS( dist_base_name, arg_1, arg_2, arg_3 ) Utility::UnitAware ## dist_base_name<arg_1,arg_2,arg_3>
%}
#define _BI_DIST_NAME_3_ARGS( dist_base_name, arg_1, arg_2, arg_3 ) Utility::UnitAware ## dist_base_name<arg_1,arg_2,arg_3>

// Get the name of a distribution that has 2 template parameters
#define _BI_DIST_NAME_2_ARGS( dist_base_name, arg_1, arg_2 ) Utility::UnitAware ## dist_base_name<arg_1,arg_2>

%{
#define BI_DIST_BASE_NAME( dist_base_name ) Utility::UnitAware ## dist_base_name
%}
// Get theunit aware base name of a distribution
#define BI_DIST_BASE_NAME( dist_base_name ) Utility::UnitAware ## dist_base_name

#define SWIG_PTR_NAME( dist_base_name ) SWIGTYPE_p_std__shared_ptrT_Utility__ ## UnitAware ## dist_base_name ## T_void_void_t_t

// Use Preprocessor Metaprogramming to get the distribution name
#define _GET_BI_DIST_NAME_MACRO( _1, _2, _3, _4, _5, _6, _7, NAME, ... ) NAME
#define BI_DIST_NAME( ... ) _GET_BI_DIST_NAME_MACRO(__VA_ARGS__, _BI_DIST_NAME_6_ARGS, _BI_DIST_NAME_5_ARGS, _BI_DIST_NAME_4_ARGS, _BI_DIST_NAME_3_ARGS, _BI_DIST_NAME_2_ARGS)(__VA_ARGS__)

//---------------------------------------------------------------------------//
// Helper macro for setting up a basic BivariateDistribution class python interface
//---------------------------------------------------------------------------//
%define %basic_bivariate_distribution_interface_setup_helper( RENAMED_DISTRIBUTION, DISTRIBUTION, PARAMS... )

%shared_ptr( BI_DIST_NAME( DISTRIBUTION, PARAMS ) )

%feature("docstring") BI_DIST_NAME( DISTRIBUTION, PARAMS )
"The RENAMED_DISTRIBUTION proxy class. This class can be evaluated and sampled.
Before sampling, make sure to initialize the Frensie Pseudo-Random Number
Generator (PyFrensie.Utility.initFrensiePrng())"

%feature("autodoc",
"evaluate(RENAMED_DISTRIBUTION self, double primary_indep_var_value, double secondary_indep_var_value ) -> double" )
BI_DIST_NAME( DISTRIBUTION, PARAMS )::evaluate;

%feature("autodoc",
"evaluateSecondaryConditionalPDF(RENAMED_DISTRIBUTION self, double primary_indep_var_value, double secondary_indep_var_value ) -> double" )
BI_DIST_NAME( DISTRIBUTION, PARAMS )::evaluateSecondaryConditionalPDF;

%feature("autodoc",
"sampleSecondaryConditional(RENAMED_DISTRIBUTION self, double primary_indep_var_value) -> double" )
BI_DIST_NAME( DISTRIBUTION, PARAMS )::sampleSecondaryConditional;

%feature("autodoc",
"sampleSecondaryConditionalAndRecordTrials(RENAMED_DISTRIBUTION self, double primary_indep_var_value, unsigned int & trials ) -> double, unsigned int

Sample from the RENAMED_DISTRIBUTION and record the number of trials
(iterations) that were required to generate the sample. The first element of
the returned tuple is the sample. The second element of the returned tuple is
the trials counter. The trials counter can be reused as such:

  PyFrensie.Utility.initFrensiePrng()

  u = PyFrensie.Utility.RENAMED_DISTRIBUTION( ... )
  trials_counter = 0

  sample,trials_counter = u.sampleSecondaryConditionalAndRecordTrials( trials_counter )
  sample,trials_counter = u.ssampleSecondaryConditionalAndRecordTrials( trials_counter )

  print trials_counter
  2")
BI_DIST_NAME( DISTRIBUTION, PARAMS )::sampleSecondaryConditionalAndRecordTrials;

%feature("autodoc",
"getUpperBoundOfSecondaryConditionalIndepVar(RENAMED_DISTRIBUTION self, double primary_indep_var_value) -> double" )
BI_DIST_NAME( DISTRIBUTION, PARAMS )::getUpperBoundOfSecondaryConditionalIndepVar;

%feature("autodoc",
"getLowerBoundOfSecondaryConditionalIndepVar(RENAMED_DISTRIBUTION self, double primary_indep_var_value) -> double" )
BI_DIST_NAME( DISTRIBUTION, PARAMS )::getLowerBoundOfSecondaryConditionalIndepVar;

// SWIG will not parse typedefs. Create some typemaps that map the
// typedefs to their true type (double)
%typemap(in) BI_DIST_NAME( DISTRIBUTION, PARAMS )::PrimaryIndepQuantity
{
  $1 = PyFloat_AsDouble($input);
}

%typemap(typecheck, precedence=90) (BI_DIST_NAME( DISTRIBUTION, PARAMS )::PrimaryIndepQuantity) {
  $1 = (PyFloat_Check($input) || PyInt_Check($input) || PyLong_Check($input)) ? 1 : 0;
}

%typemap(out) BI_DIST_NAME( DISTRIBUTION, PARAMS )::PrimaryIndepQuantity
{
  $result = PyFloat_FromDouble($1);
}

%typemap(in) BI_DIST_NAME( DISTRIBUTION, PARAMS )::SecondaryIndepQuantity
{
  $1 = PyFloat_AsDouble($input);
}

%typemap(typecheck, precedence=90) (BI_DIST_NAME( DISTRIBUTION, PARAMS )::SecondaryIndepQuantity) {
  $1 = (PyFloat_Check($input) || PyInt_Check($input) || PyLong_Check($input)) ? 1 : 0;
}

%typemap(out) BI_DIST_NAME( DISTRIBUTION, PARAMS )::SecondaryIndepQuantity
{
  $result = PyFloat_FromDouble($1);
}

%typemap(out) BI_DIST_NAME( DISTRIBUTION, PARAMS )::DepQuantity
{
  $result = PyFloat_FromDouble($1);
}

%typemap(out) BI_DIST_NAME( DISTRIBUTION, PARAMS )::InverseSecondaryIndepQuantity
{
  $result = PyFloat_FromDouble($1);
}

%template(RENAMED_DISTRIBUTION) BI_DIST_NAME( DISTRIBUTION, PARAMS );

%enddef

//---------------------------------------------------------------------------//
// Helper macro for setting up a basic TabularBivariateDistribution class py. int.
//---------------------------------------------------------------------------//
%define %standard_tab_bivariate_distribution_interface_setup_helper( RENAMED_DISTRIBUTION, DISTRIBUTION, PARAMS... )

%feature("autodoc",
"evaluateSecondaryConditionalCDF(RENAMED_DISTRIBUTION self, double primary_indep_var_value, double secondary_indep_var_value ) -> double" )
BI_DIST_NAME( DISTRIBUTION, PARAMS )::evaluateSecondaryConditionalCDF;

%feature("autodoc",
"sampleSecondaryConditionalAndRecordBinIndices(RENAMED_DISTRIBUTION self, double primary_indep_var_value) -> double, size_t, size_t

Sample from the RENAMED_DISTRIBUTION and record the bin index corresponding to
the sample. The first element of the returned tuple is the sample. The second
and third elements of the returned tuple are the primary and secondary bin
indices. This method can be called as follows:

  PyFrensie.Utility.initFrensiePrng()

  u = PyFrensie.Utility.RENAMED_DISTRIBUTION( ... )
  primary_indep_var_value = 1.0

  sample,primary_bin_index,secondary_bin_index = u.sampleSecondaryConditionalAndRecordBinIndices( primary_indep_var_value)")
BI_DIST_NAME( DISTRIBUTION, PARAMS )::sampleSecondaryConditionalAndRecordBinIndices;

%feature("autodoc",
"sampleSecondaryConditionalWithRandomNumber(RENAMED_DISTRIBUTION self, double primary_indep_var_value, const double random_number) -> double

Sample from the RENAMED_DISTRIBUTION using the supplied random number instead
of using the hidden Utility::RandomNumberGenerator.")
BI_DIST_NAME( DISTRIBUTION, PARAMS )::sampleSecondaryConditionalWithRandomNumber;

%feature("autodoc",
"sampleSecondaryConditionalInSubrange(RENAMED_DISTRIBUTION self, double primary_indep_var_value, const double max_secondary_indep_var_value ) -> double

Sample from the RENAMED_DISTRIBUTION in the subrange
[self.getLowerBoundOfSecondaryConditionalIndepVar( primary_indep_var_value ),max_secondary_indep_var_value]" )
BI_DIST_NAME( DISTRIBUTION, PARAMS )::sampleSecondaryConditionalInSubrange;

%feature("autodoc",
"sampleSecondaryConditionalWithRandomNumberInSubrange(RENAMED_DISTRIBUTION self, double primary_indep_var_value, const double random_number, const double max_secondary_indep_var_value ) -> double

Sample from the RENAMED_DISTRIBUTION using the supplied random number in the
subrange [self.getLowerBoundOfSecondaryConditionalIndepVar( primary_indep_var_value ),max_indep_var]" )
BI_DIST_NAME( DISTRIBUTION, PARAMS )::sampleSecondaryConditionalWithRandomNumberInSubrange;

// // Create unitless constructor
// %extend BI_DIST_BASE_NAME( DISTRIBUTION )
// {
  // // Constructor
  // BI_DIST_BASE_NAME( DISTRIBUTION )(
  //   const std::vector<double>& raw_prim_grid,
  //   const std::vector<std::shared_ptr<Utility::UnitAwareTabularUnivariateDistribution<void,void> > >& vec )
  // {
  //   std::vector<double> prim_grid(raw_prim_grid);
  //   std::vector<std::shared_ptr<const Utility::UnitAwareTabularUnivariateDistribution<void,void> > > sec_dists(vec.size());

  //   for ( unsigned i = 0; i < vec.size(); ++i )
  //   {
  //     sec_dists[i] = vec[i];
  //   }

  //   return new BI_DIST_NAME( DISTRIBUTION, PARAMS )( prim_grid, sec_dists );
  // }
// }

%basic_bivariate_distribution_interface_setup_helper( RENAMED_DISTRIBUTION, DISTRIBUTION, PARAMS )

%enddef

//---------------------------------------------------------------------------//
// Helper macros for extending a BivariateDistribution class python interface
//---------------------------------------------------------------------------//
%define %extend_tab_bivariate_distribution_interface_helper( RENAMED_DISTRIBUTION, DISTRIBUTION, BASE_DISTRIBUTION, PARAMS... )

// Typemap for std::vector<PrimaryIndepQuantity
%typemap(in) const std::vector<BI_DIST_NAME( DISTRIBUTION, PARAMS )::PrimaryIndepQuantity,std::allocator< BI_DIST_NAME( DISTRIBUTION, PARAMS )::PrimaryIndepQuantity > > & (std::vector<double> temp)
{
  temp = PyFrensie::convertFromPython<std::vector<double> >( $input );

  $1 = &temp;
}

// Typecheck for std::vector<PrimaryIndepQuantity>
%typemap(typecheck, precedence=1050) (const std::vector< BI_DIST_NAME( DISTRIBUTION, PARAMS )::PrimaryIndepQuantity,std::allocator< BI_DIST_NAME( DISTRIBUTION, PARAMS )::PrimaryIndepQuantity > > &) {
  $1 = (PyArray_Check($input) || PySequence_Check($input)) ? 1 : 0;
}

// Typemap for std::vector< std::vector< SecondaryIndepQuantity > >
%typemap(in) const std::vector< std::vector< BI_DIST_NAME( DISTRIBUTION, PARAMS )::SecondaryIndepQuantity,std::allocator< BI_DIST_NAME( DISTRIBUTION, PARAMS )::SecondaryIndepQuantity > >,std::allocator< std::vector< BI_DIST_NAME( DISTRIBUTION, PARAMS )::SecondaryIndepQuantity,std::allocator< BI_DIST_NAME( DISTRIBUTION, PARAMS )::SecondaryIndepQuantity > > > > & (std::vector<std::vector<double> > temp)
{
  temp = PyFrensie::convertFromPython<std::vector<std::vector<double> > >( $input );

  $1 = &temp;
}

// Typecheck for std::vector< std::vector< SecondaryIndepQuantity > >
%typemap(typecheck, precedence=1050) (const std::vector< std::vector< BI_DIST_NAME( DISTRIBUTION, PARAMS )::SecondaryIndepQuantity,std::allocator< BI_DIST_NAME( DISTRIBUTION, PARAMS )::SecondaryIndepQuantity > >,std::allocator< std::vector< BI_DIST_NAME( DISTRIBUTION, PARAMS )::SecondaryIndepQuantity,std::allocator< BI_DIST_NAME( DISTRIBUTION, PARAMS )::SecondaryIndepQuantity > > > > &) {
  $1 = (PyArray_Check($input) || PySequence_Check($input)) ? 1 : 0;
}

// Typemap for std::vector< std::vector< DepQuantity > >
%typemap(in) const std::vector< std::vector< BI_DIST_NAME( DISTRIBUTION, PARAMS )::DepQuantity,std::allocator< BI_DIST_NAME( DISTRIBUTION, PARAMS )::DepQuantity > >,std::allocator< std::vector< BI_DIST_NAME( DISTRIBUTION, PARAMS )::DepQuantity,std::allocator< BI_DIST_NAME( DISTRIBUTION, PARAMS )::DepQuantity > > > > & (std::vector<std::vector<double> > temp)
{
  temp = PyFrensie::convertFromPython<std::vector<std::vector<double> > >( $input );

  $1 = &temp;
}

// Typecheck for std::vector< std::vector< DepQuantity > >
%typemap(typecheck, precedence=1050) (const std::vector< std::vector< BI_DIST_NAME( DISTRIBUTION, PARAMS )::DepQuantity,std::allocator< BI_DIST_NAME( DISTRIBUTION, PARAMS )::DepQuantity > >,std::allocator< std::vector< BI_DIST_NAME( DISTRIBUTION, PARAMS )::DepQuantity,std::allocator< BI_DIST_NAME( DISTRIBUTION, PARAMS )::DepQuantity > > > > &) {
  $1 = (PyArray_Check($input) || PySequence_Check($input)) ? 1 : 0;
}

%extend BI_DIST_NAME( DISTRIBUTION, PARAMS )
{
  // Evaluate the distribution
  DepQuantity BI_DIST_NAME( DISTRIBUTION, PARAMS )::evaluate(
            const PrimaryIndepQuantity primary_indep_var_value,
            const SecondaryIndepQuantity secondary_indep_var_value ) const
  {
    return $self->evaluate( primary_indep_var_value, secondary_indep_var_value );
  }

  // Evaluate the secondary conditional PDF
  InverseSecondaryIndepQuantity BI_DIST_NAME( DISTRIBUTION, PARAMS )::evaluateSecondaryConditionalPDF(
            const PrimaryIndepQuantity primary_indep_var_value,
            const SecondaryIndepQuantity secondary_indep_var_value ) const
  {
    return $self->evaluateSecondaryConditionalPDF( primary_indep_var_value,
                                                   secondary_indep_var_value );
  }

  // Return the upper bound of the distribution primary independent variable
  PrimaryIndepQuantity BI_DIST_NAME( DISTRIBUTION, PARAMS )::getUpperBoundOfPrimaryIndepVar() const
  {
    return $self->getUpperBoundOfPrimaryIndepVar();
  }

  // Return the lower bound of the distribution primary independent variable
  PrimaryIndepQuantity BI_DIST_NAME( DISTRIBUTION, PARAMS )::getLowerBoundOfPrimaryIndepVar() const
  {
    return $self->getLowerBoundOfPrimaryIndepVar();
  }

  // Return the upper bound of the secondary conditional distribution
  SecondaryIndepQuantity BI_DIST_NAME( DISTRIBUTION, PARAMS )::getUpperBoundOfSecondaryConditionalIndepVar(
                const PrimaryIndepQuantity primary_indep_var_value ) const
  {
    return $self->getUpperBoundOfSecondaryConditionalIndepVar( primary_indep_var_value );
  }

  // Return the lower bound of the secondary conditional distribution
  SecondaryIndepQuantity BI_DIST_NAME( DISTRIBUTION, PARAMS )::getLowerBoundOfSecondaryConditionalIndepVar(
                const PrimaryIndepQuantity primary_indep_var_value ) const
  {
    return $self->getLowerBoundOfSecondaryConditionalIndepVar( primary_indep_var_value );
  }

  // Test if the distribution is tabular in the primary dimension
  bool BI_DIST_NAME( DISTRIBUTION, PARAMS )::isPrimaryDimensionTabular() const
  {
    return $self->isPrimaryDimensionTabular();
  }

  // Test if the distribution is continuous in the primary dimension
  bool BI_DIST_NAME( DISTRIBUTION, PARAMS )::isPrimaryDimensionContinuous() const
  {
    return $self->isPrimaryDimensionContinuous();
  }

  // // Test if the distribution has the same primary bounds
  // bool BI_DIST_NAME( DISTRIBUTION, PARAMS )::hasSamePrimaryBounds( UnitAwareBasicBivariateDistribution<void,void,void>& distribution ) const
  // {
  //   return $self->hasSamePrimaryBounds( distribution );
  // }

};

%enddef

//---------------------------------------------------------------------------//
// Macro for setting up a basic Bivariate Distribution class python interface
//---------------------------------------------------------------------------//
%define %basic_bivariate_distribution_interface_setup( DISTRIBUTION )

%basic_bivariate_distribution_interface_setup_helper( DISTRIBUTION, DISTRIBUTION, void, void, void )

%enddef

//---------------------------------------------------------------------------//
// Macro for setting up a standard Bivariate Distribution class python interface
//---------------------------------------------------------------------------//
%define %standard_bivariate_distribution_interface_setup( DISTRIBUTION )

// Do the basic setup for this distribution
%basic_bivariate_distribution_interface_setup_helper( DISTRIBUTION, DISTRIBUTION, void ,void, void )

%enddef

//---------------------------------------------------------------------------//
// Macro for setting up an advanced Bivariate Distribution class python interface
//---------------------------------------------------------------------------//
%define %advanced_bivariate_distribution_interface_setup( RENAMED_DISTRIBUTION, DISTRIBUTION, INTERP, GRID, BASE_DISTRIBUTION )

// Add a typedef for the renamed distribution so that the extended methods
// can be compiled
%inline %{
typedef BI_DIST_NAME( DISTRIBUTION, INTERP, GRID, BASE_DISTRIBUTION, void, void, void ) RENAMED_DISTRIBUTION;
%}

// Do the basic setup for this distribution
%basic_bivariate_distribution_interface_setup_helper( RENAMED_DISTRIBUTION, DISTRIBUTION, INTERP, GRID, BASE_DISTRIBUTION, void, void, void )

%enddef

//---------------------------------------------------------------------------//
// Macro for setting up a basic Tabular Bivariate Distribution class py. int.
//---------------------------------------------------------------------------//
%define %basic_tab_bivariate_distribution_interface_setup( DISTRIBUTION, PARAMETER )

// Add a typedef for the renamed distribution so that the extended methods
// can be compiled
%inline %{
typedef BI_DIST_NAME( DISTRIBUTION, void, void, void, PARAMETER ) DISTRIBUTION;
%}

// Do the basic tabular setup for this distribution
%basic_bivariate_distribution_interface_setup_helper( DISTRIBUTION, DISTRIBUTION, void , void, void, PARAMETER )

%enddef

//---------------------------------------------------------------------------//
// Macro for setting up a standard Tabular Bivariate Distribution class python interface
//---------------------------------------------------------------------------//
%define %standard_tab_bivariate_distribution_interface_setup( DISTRIBUTION )

%standard_tab_bivariate_distribution_interface_setup_helper( DISTRIBUTION, DISTRIBUTION, void, void, void )

%enddef

//---------------------------------------------------------------------------//
// Macro for setting up an advanced Tabular Bivariate Distribution class python interface
//---------------------------------------------------------------------------//
%define %advanced_tab_bivariate_distribution_interface_setup( RENAMED_DISTRIBUTION, DISTRIBUTION, BASE_DISTRIBUTION, INTERP, GRID )

// Add a typedef for the renamed distribution so that the extended methods
// can be compiled
%inline %{
typedef BI_DIST_NAME( DISTRIBUTION, INTERP, GRID, void, void, void ) RENAMED_DISTRIBUTION;
%}

%extend BI_DIST_NAME( DISTRIBUTION, INTERP, GRID, void, void, void )
{
  // Constructor
  BI_DIST_NAME( DISTRIBUTION, INTERP, GRID, void, void, void )(
    const std::vector<double>& prim_grid,
    PyObject* raw_sec_dists,
    const double fuzzy_boundary_tol = 1e-3,
    const double evaluate_relative_error_tol = 1e-7,
    const double evaluate_error_tol = 1e-16 )
  {
    std::vector<std::shared_ptr<const BI_DIST_NAME(BASE_DISTRIBUTION, void, void ) > > sec_dists;

    bool conversion_success;

    CONVERT_PYOBJECT_TO_VECTOR_OF_BASE_SHARED_PTR( raw_sec_dists, sec_dists, SWIG_PTR_NAME( BASE_DISTRIBUTION ), conversion_success );

    if( conversion_success )
    {
      return new BI_DIST_NAME( DISTRIBUTION, INTERP, GRID, void, void, void )(
            prim_grid,
            sec_dists,
            fuzzy_boundary_tol,
            evaluate_relative_error_tol,
            evaluate_error_tol );
    }
    // SWIG will check for a NULL return type and throw an exception
    else
      return NULL;
  }
};

%extend_tab_bivariate_distribution_interface_helper( RENAMED_DISTRIBUTION, DISTRIBUTION, BASE_DISTRIBUTION, INTERP, GRID, void, void, void )

// Do the basic tabular setup for this distribution
%standard_tab_bivariate_distribution_interface_setup_helper( RENAMED_DISTRIBUTION, DISTRIBUTION, INTERP, GRID, void, void, void )

%enddef

//---------------------------------------------------------------------------//
// Macro for setting up an advanced Histogram Tabular Bivariate Distribution class python interface
//---------------------------------------------------------------------------//
%define %advanced_histogram_bivariate_distribution_interface_setup( RENAMED_DISTRIBUTION, DISTRIBUTION, BASE_DISTRIBUTION )

// Add a typedef for the renamed distribution so that the extended methods
// can be compiled
%inline %{
typedef BI_DIST_NAME( DISTRIBUTION, void, void, void ) RENAMED_DISTRIBUTION;
%}

%extend BI_DIST_NAME( DISTRIBUTION, void, void, void )
{
  // Constructor
  BI_DIST_NAME( DISTRIBUTION, void, void, void )(
    const std::vector<double>& prim_grid,
    PyObject* raw_sec_dists )
  {
    std::vector<std::shared_ptr<const BI_DIST_NAME(BASE_DISTRIBUTION, void, void ) > > sec_dists;

    bool conversion_success;

    CONVERT_PYOBJECT_TO_VECTOR_OF_BASE_SHARED_PTR( raw_sec_dists, sec_dists, SWIG_PTR_NAME( BASE_DISTRIBUTION ), conversion_success );

    if( conversion_success )
    {
      return new BI_DIST_NAME( DISTRIBUTION, void, void, void )( prim_grid, sec_dists );
    }
    // SWIG will check for a NULL return type and throw an exception
    else
      return NULL;
  }

  // Return a random sample from the secondary conditional PDF
  SecondaryIndepQuantity BI_DIST_NAME( DISTRIBUTION, void, void, void )::sampleSecondaryConditional(
                const PrimaryIndepQuantity primary_indep_var_value ) const
  {
    return $self->sampleSecondaryConditional( primary_indep_var_value );
  }

  // Return a random sample and record the number of trials
  SecondaryIndepQuantity BI_DIST_NAME( DISTRIBUTION, void, void, void )::sampleSecondaryConditionalAndRecordTrials(
                            const PrimaryIndepQuantity primary_indep_var_value,
                            DistributionTraits::Counter& trials ) const
  {
    return $self->sampleSecondaryConditionalAndRecordTrials( primary_indep_var_value, trials);
  }
};

%extend_tab_bivariate_distribution_interface_helper( RENAMED_DISTRIBUTION, DISTRIBUTION, BASE_DISTRIBUTION, void, void, void )

// Do the basic tabular setup for this distribution
%standard_tab_bivariate_distribution_interface_setup_helper( RENAMED_DISTRIBUTION, DISTRIBUTION, void, void, void )

%enddef

//---------------------------------------------------------------------------//
// end Utility_BivariateDistributionHelpers.i
//---------------------------------------------------------------------------//
