//---------------------------------------------------------------------------//
//!
//! \file   tstList.cpp
//! \author Alex Robinson
//! \brief  List helper function unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>
#include <sstream>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>

// FRENSIE Includes
#include "Utility_List.hpp"
#include "Utility_Tuple.hpp"
#include "Utility_UnitTestHarnessExtensions.hpp"

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that a list can be converted to a string
TEUCHOS_UNIT_TEST( List, toString )
{
  TEST_EQUALITY_CONST( Utility::toString( std::list<short>( {-1, 2} ) ),
                       "{-1, 2}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<unsigned short>( {0, 10, 100} ) ),
                       "{0, 10, 100}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<int>( {-11111, 0, 11111, 22222} ) ),
                       "{-11111, 0, 11111, 22222}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<unsigned int>( {0, 10, 100, 1000, 10000} ) ),
                       "{0, 10, 100, 1000, 10000}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<long>( {-11111, 0, 11111, 22222} ) ),
                       "{-11111, 0, 11111, 22222}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<unsigned long>( {0, 10, 100, 1000, 10000} ) ),
                       "{0, 10, 100, 1000, 10000}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<long long>( {-10000000000, 0, 10000000000} ) ),
                       "{-10000000000, 0, 10000000000}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<unsigned long long>( {0, 1000000000, 1000000000000} ) ),
                       "{0, 1000000000, 1000000000000}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<float>( {-1.0f, 0.0f, 1.0f} ) ),
                       "{-1.000000000e+00, 0.000000000e+00, 1.000000000e+00}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<double>( {-1.0, 0.0, 1.0} ) ),
                       "{-1.000000000000000000e+00, 0.000000000000000000e+00, 1.000000000000000000e+00}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<char>( {'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g' } ) ),
                       "{T, e, s, t,  , s, t, r, i, n, g}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<std::string>( {"Test", "string"} ) ),
                       "{Test, string}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<std::pair<int, int> >( {std::pair<int, int>({0, 1}), std::pair<int,int>({-1, 2})} ) ),
                       "{{0, 1}, {-1, 2}}" );
  TEST_EQUALITY_CONST( Utility::toString( std::list<std::tuple<unsigned, double, long> >( {std::make_tuple( 0u, 1.0, -100000 ), std::make_tuple( 1u, -1.0, 100001 )} ) ),
                       "{{0, 1.000000000000000000e+00, -100000}, {1, -1.000000000000000000e+00, 100001}}" );
}

//---------------------------------------------------------------------------//
// Check that a list can be placed in a stream
TEUCHOS_UNIT_TEST( List, toStream )
{
  std::ostringstream oss;

  Utility::toStream( oss, std::list<short>( {-1, 2} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{-1, 2}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<unsigned short>( {0, 10, 100} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<int>( {-11111, 0, 11111, 22222} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{-11111, 0, 11111, 22222}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<unsigned int>( {0, 10, 100, 1000, 10000} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100, 1000, 10000}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<long>( {-11111, 0, 11111, 22222} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{-11111, 0, 11111, 22222}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<unsigned long>( {0, 10, 100, 1000, 10000} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100, 1000, 10000}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<long long>( {-10000000000, 0, 10000000000} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{-10000000000, 0, 10000000000}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<unsigned long long>( {0, 1000000000, 1000000000000} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 1000000000, 1000000000000}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<float>( {-1.0f, 0.0f, 1.0f} ) );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{-1.000000000e+00, 0.000000000e+00, 1.000000000e+00}" );
  
  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<double>( {-1.0, 0.0, 1.0} ) );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{-1.000000000000000000e+00, 0.000000000000000000e+00, 1.000000000000000000e+00}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<char>( {'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g' } ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{T, e, s, t,  , s, t, r, i, n, g}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<std::string>( {"Test", "string"} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{Test, string}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<std::pair<int, int>>( {std::pair<int, int>({0, 1}), std::pair<int,int>({-1, 2})} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{{0, 1}, {-1, 2}}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::list<std::tuple<unsigned, double, long>>( {std::make_tuple( 0u, 1.0, -100000 ), std::make_tuple( 1u, -1.0, 100001 )} ) );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{{0, 1.000000000000000000e+00, -100000}, {1, -1.000000000000000000e+00, 100001}}" );
}

//---------------------------------------------------------------------------//
// Check that a forward_list can be converted to a string
TEUCHOS_UNIT_TEST( ForwardList, toString )
{
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<short>( {-1, 2} ) ),
                       "{-1, 2}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<unsigned short>( {0, 10, 100} ) ),
                       "{0, 10, 100}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<int>( {-11111, 0, 11111, 22222} ) ),
                       "{-11111, 0, 11111, 22222}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<unsigned int>( {0, 10, 100, 1000, 10000} ) ),
                       "{0, 10, 100, 1000, 10000}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<long>( {-11111, 0, 11111, 22222} ) ),
                       "{-11111, 0, 11111, 22222}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<unsigned long>( {0, 10, 100, 1000, 10000} ) ),
                       "{0, 10, 100, 1000, 10000}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<long long>( {-10000000000, 0, 10000000000} ) ),
                       "{-10000000000, 0, 10000000000}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<unsigned long long>( {0, 1000000000, 1000000000000} ) ),
                       "{0, 1000000000, 1000000000000}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<float>( {-1.0f, 0.0f, 1.0f} ) ),
                       "{-1.000000000e+00, 0.000000000e+00, 1.000000000e+00}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<double>( {-1.0, 0.0, 1.0} ) ),
                       "{-1.000000000000000000e+00, 0.000000000000000000e+00, 1.000000000000000000e+00}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<char>( {'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g' } ) ),
                       "{T, e, s, t,  , s, t, r, i, n, g}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<std::string>( {"Test", "string"} ) ),
                       "{Test, string}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<std::pair<int, int> >( {std::pair<int, int>({0, 1}), std::pair<int,int>({-1, 2})} ) ),
                       "{{0, 1}, {-1, 2}}" );
  TEST_EQUALITY_CONST( Utility::toString( std::forward_list<std::tuple<unsigned, double, long> >( {std::make_tuple( 0u, 1.0, -100000 ), std::make_tuple( 1u, -1.0, 100001 )} ) ),
                       "{{0, 1.000000000000000000e+00, -100000}, {1, -1.000000000000000000e+00, 100001}}" );
}

//---------------------------------------------------------------------------//
// Check that a forward_list can be placed in a stream
TEUCHOS_UNIT_TEST( ForwardList, toStream )
{
  std::ostringstream oss;

  Utility::toStream( oss, std::forward_list<short>( {-1, 2} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{-1, 2}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<unsigned short>( {0, 10, 100} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<int>( {-11111, 0, 11111, 22222} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{-11111, 0, 11111, 22222}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<unsigned int>( {0, 10, 100, 1000, 10000} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100, 1000, 10000}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<long>( {-11111, 0, 11111, 22222} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{-11111, 0, 11111, 22222}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<unsigned long>( {0, 10, 100, 1000, 10000} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100, 1000, 10000}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<long long>( {-10000000000, 0, 10000000000} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{-10000000000, 0, 10000000000}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<unsigned long long>( {0, 1000000000, 1000000000000} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 1000000000, 1000000000000}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<float>( {-1.0f, 0.0f, 1.0f} ) );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{-1.000000000e+00, 0.000000000e+00, 1.000000000e+00}" );
  
  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<double>( {-1.0, 0.0, 1.0} ) );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{-1.000000000000000000e+00, 0.000000000000000000e+00, 1.000000000000000000e+00}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<char>( {'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g' } ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{T, e, s, t,  , s, t, r, i, n, g}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<std::string>( {"Test", "string"} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{Test, string}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<std::pair<int, int>>( {std::pair<int, int>({0, 1}), std::pair<int,int>({-1, 2})} ) );
  
  TEST_EQUALITY_CONST( oss.str(), "{{0, 1}, {-1, 2}}" );

  oss.str( "" );
  oss.clear();

  Utility::toStream( oss, std::forward_list<std::tuple<unsigned, double, long>>( {std::make_tuple( 0u, 1.0, -100000 ), std::make_tuple( 1u, -1.0, 100001 )} ) );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{{0, 1.000000000000000000e+00, -100000}, {1, -1.000000000000000000e+00, 100001}}" );
}

//---------------------------------------------------------------------------//
// Check that a list can be created from a string
TEUCHOS_UNIT_TEST( List, fromString )
{
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<short> >( "{-1, 2}" )),
                           std::list<short>({-1, 2}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<unsigned short> >( "{0, 10, 100}" )),
                           std::list<unsigned short>({0, 10, 100}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<int> >( "{-11111, 0, 11111, 22222}" )),
                           std::list<int>({-11111, 0, 11111, 22222}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<unsigned int> >( "{0, 10, 100, 1000, 10000}" )),
                           std::list<unsigned int>({0, 10, 100, 1000, 10000}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<long> >( "{-11111, 0, 11111, 22222}" )),
                           std::list<long>({-11111, 0, 11111, 22222}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<unsigned long> >( "{0, 10, 100, 1000, 10000}" )),
                           std::list<unsigned long>({0, 10, 100, 1000, 10000}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<long long> >( "{-1000000000, 0, 1000000000}" )),
                           std::list<long long>({-1000000000, 0, 1000000000}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<unsigned long long> >( "{0, 1000000000, 10000000000}" )),
                           std::list<unsigned long long>({0, 1000000000, 10000000000}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<float> >( "{-1, 0.0, 1.000000000e+00}" )),
                           std::list<float>({-1.0f, 0.0f, 1.0f}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<double> >( "{-1, 0.0, 1.000000000000000000e+00}" )),
                           std::list<double>({-1.0, 0.0, 1.0}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<char> >( "{T, e, s, t,  , s, t, r, i, n, g}" )),
                           std::list<char>({'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g'}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<std::string> >( "{Test, string}" )),
                           std::list<std::string>({"Test","string"}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<std::pair<int, int> > >( "{{0, 1}, {-1, 2}}" )),
                           (std::list<std::pair<int,int> >({std::make_pair(0, 1), std::make_pair(-1, 2)})) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::list<std::tuple<unsigned, double, long> > >( "{{0, 1.0, -100000}, {1, -1.00, 100001}}" )),
                           (std::list<std::tuple<unsigned,double,long> >({std::make_tuple(0u, 1.0, -100000l), std::make_tuple(1u, -1.0, 100001l)})) );
}

//---------------------------------------------------------------------------//
// Check that a list can be extracted from a stream
TEUCHOS_UNIT_TEST( List, fromStream )
{
  // Extract list of short
  std::istringstream iss( "{-1, 2}" );

  {
    std::list<short> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<short>({-1, 2}) );

    iss.str( "{-1, 2i, 2}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<short>({-1, 0, 1, 2}) );
  }

  // Extract list of unsigned short
  iss.str( "{0, 10, 100}" );
  iss.clear();

  {
    std::list<unsigned short> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned short>({0, 10, 100}) );

    iss.str( "{0, 3i, 100}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned short>({0, 25, 50, 75, 100}) );
  }

  // Extract list of int
  iss.str( "{-11111, 0, 11111, 22222}" );
  iss.clear();

  {
    std::list<int> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<int>({-11111, 0, 11111, 22222}) );

    iss.str( "{-11111, 2i, 22222}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<int>({-11111, 0, 11111, 22222}) );
  }

  // Extract list of unsigned int
  iss.str( "{0, 10, 100, 1000, 10000}" );
  iss.clear();

  {
    std::list<unsigned> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned>({0, 10, 100, 1000, 10000}) );

    iss.str( "{0, 3i, 10000}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned>({0, 2500, 5000, 7500, 10000}) );
  }

  // Extract list of long
  iss.str( "{-11111, 0, 11111, 22222}" );
  iss.clear();

  {
    std::list<long> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<long>({-11111, 0, 11111, 22222}) );

    iss.str( "{-11111, 2i, 22222}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<long>({-11111, 0, 11111, 22222}) );
  }

  // Extract list of unsigned long
  iss.str( "{0, 10, 100, 1000, 10000}" );
  iss.clear();

  {
    std::list<unsigned long> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned long>({0, 10, 100, 1000, 10000}) );

    iss.str( "{0, 3i, 10000}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned long>({0, 2500, 5000, 7500, 10000}) );
  }

  // Extract list of long long
  iss.str( "{-1000000000, 0, 1000000000}" );
  iss.clear();

  {
    std::list<long long> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<long long>({-1000000000, 0, 1000000000}) );

    iss.str( "{-1000000000, 1i, 1000000000}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<long long>({-1000000000, 0, 1000000000}) );
  }

  // Extract list of unsigned long long
  iss.str( "{0, 1000000000, 10000000000}" );
  iss.clear();

  {
    std::list<unsigned long long> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned long long>({0, 1000000000, 10000000000}) );

    iss.str( "{0, 1i, 10000000000}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned long long>({0, 5000000000, 10000000000}) );
  }

  // Extract list of float
  iss.str( "{-1, 0.0, 1.000000000e+00}" );
  iss.clear();

  {
    std::list<float> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<float>({-1.0f, 0.0f, 1.0f}) );

    iss.str( "{-1, 3i, 1.000000000e+00}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<float>({-1.0f, -0.5f, 0.0f, 0.5f, 1.0f}) );

    iss.str( "{1e-3, 2l, 1.0}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::list<float>({1e-3f, 1e-2f, 1e-1f, 1.0f}),
                                      1e-7 );

    iss.str( "{1.0, 1l, 100, 3i, 200}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<float>({1.0f, 10.0f, 100.0f, 125.0f, 150.0f, 175.0f, 200.0f}) );

    iss.str( "{0.0, pi/2, 3Pi / 4, PI, 2*pi}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<float>({0.0f, (float)Utility::PhysicalConstants::pi/2, 3*(float)Utility::PhysicalConstants::pi/4, (float)Utility::PhysicalConstants::pi, 2*(float)Utility::PhysicalConstants::pi}) );

    iss.str( "{-pi, 3i, 0}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::list<float>({-(float)Utility::PhysicalConstants::pi, -3*(float)Utility::PhysicalConstants::pi/4, -(float)Utility::PhysicalConstants::pi/2, -(float)Utility::PhysicalConstants::pi/4, 0.0}),
                                      1e-7 );

    iss.str( "{-inf, 0, Infinity}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<float>({-std::numeric_limits<float>::infinity(), 0.0f, std::numeric_limits<float>::infinity()}) );
  }

  // Extract list of double
  iss.str( "{-1, 0.0, 1.000000000000000000e+00}" );
  iss.clear();

  {
    std::list<double> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<double>({-1.0, 0.0, 1.0}) );

    iss.str( "{-1, 3i, 1.000000000e+00}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<double>({-1.0, -0.5, 0.0, 0.5, 1.0}) );

    iss.str( "{1e-3, 2l, 1.0}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::list<double>({1e-3, 1e-2, 1e-1, 1.0}),
                                      1e-9 );

    iss.str( "{1.0, 1l, 100, 3i, 200}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::list<double>({1.0, 10.0, 100.0, 125.0, 150.0, 175.0, 200.0}),
                                      1e-15 );

    iss.str( "{0.0, pi/2, 3Pi / 4, PI, 2*pi}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<double>({0.0, Utility::PhysicalConstants::pi/2, 3*Utility::PhysicalConstants::pi/4, Utility::PhysicalConstants::pi, 2*Utility::PhysicalConstants::pi}) );

    iss.str( "{-pi, 3i, 0}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<double>({-Utility::PhysicalConstants::pi, -3*Utility::PhysicalConstants::pi/4, -Utility::PhysicalConstants::pi/2, -Utility::PhysicalConstants::pi/4, 0.0}) );

    iss.str( "{-inf, 0, Infinity}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<double>({-std::numeric_limits<float>::infinity(), 0.0, std::numeric_limits<float>::infinity()}) );
  }

  // Extract list of char
  iss.str( "{T, e, s, t,  , s, t, r, i, n, g}" );
  iss.clear();

  {
    std::list<char> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<char>({'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g'}) );
  }

  // Extract list of string
  iss.str( "{Test, string}" );
  iss.clear();

  {
    std::list<std::string> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::list<std::string>({"Test","string"}) );
  }

  // Extract list of pair
  iss.str( "{{0, 1}, {-1, 2}}" );
  iss.clear();

  {
    std::list<std::pair<int,int> > test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, (std::list<std::pair<int,int> >({std::make_pair(0, 1), std::make_pair(-1, 2)})) );
  }

  // Extract list of tuple
  iss.str( "{{0, 1.0, -100000}, {1, -1.00, 100001}}" );
  iss.clear();

  {
    std::list<std::tuple<unsigned,double,long> > test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, (std::list<std::tuple<unsigned,double,long> >({std::make_tuple(0u, 1.0, -100000l), std::make_tuple(1u, -1.0, 100001l)})) );
  }

  // Extract multiple lists of pairs
  iss.str( "{{0, 1}, {-1, 2}}, {{1, 0}, {2, -1}}" );
  iss.clear();

  {
    std::list<std::pair<long,long> > test_list;

    Utility::fromStream( iss, test_list, "," );

    TEST_COMPARE_CONTAINERS( test_list, (std::list<std::pair<long,long> >({std::make_pair(0l, 1l), std::make_pair(-1l, 2l)})) );

    Utility::moveInputStreamToNextElement( iss, ',', '}' );

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, (std::list<std::pair<long,long> >({std::make_pair(1l, 0l), std::make_pair(2l, -1l)})) );
  }

  // Extract list of lists
  iss.str( "{{-1, 3i, 1}, {1e-3, 2l, 1.0}, {-Inf, Inf}, {-pi/2, 1i, pi/2}}" );
  iss.clear();

  {
    std::list<std::list<double> > test_list;
    
    Utility::fromStream( iss, test_list );
    
    TEST_EQUALITY_CONST( test_list.size(), 4 );

    std::list<std::list<double> >::const_iterator
      test_list_it = test_list.begin();
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-1.0, -0.5, 0.0, 0.5, 1.0})) );
    
    ++test_list_it;
    
    TEST_COMPARE_FLOATING_CONTAINERS( *test_list_it,
                                      (std::list<double>({1e-3, 1e-2, 1e-1, 1.0})),
                                      1e-9 );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()})) );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-Utility::PhysicalConstants::pi/2, 0.0, Utility::PhysicalConstants::pi/2})) );
  }
}

//---------------------------------------------------------------------------//
// Check that a forward_list can be created from a string
TEUCHOS_UNIT_TEST( ForwardList, fromString )
{
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<short> >( "{-1, 2}" )),
                           std::forward_list<short>({-1, 2}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<unsigned short> >( "{0, 10, 100}" )),
                           std::forward_list<unsigned short>({0, 10, 100}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<int> >( "{-11111, 0, 11111, 22222}" )),
                           std::forward_list<int>({-11111, 0, 11111, 22222}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<unsigned int> >( "{0, 10, 100, 1000, 10000}" )),
                           std::forward_list<unsigned int>({0, 10, 100, 1000, 10000}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<long> >( "{-11111, 0, 11111, 22222}" )),
                           std::forward_list<long>({-11111, 0, 11111, 22222}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<unsigned long> >( "{0, 10, 100, 1000, 10000}" )),
                           std::forward_list<unsigned long>({0, 10, 100, 1000, 10000}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<long long> >( "{-1000000000, 0, 1000000000}" )),
                           std::forward_list<long long>({-1000000000, 0, 1000000000}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<unsigned long long> >( "{0, 1000000000, 10000000000}" )),
                           std::forward_list<unsigned long long>({0, 1000000000, 10000000000}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<float> >( "{-1, 0.0, 1.000000000e+00}" )),
                           std::forward_list<float>({-1.0f, 0.0f, 1.0f}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<double> >( "{-1, 0.0, 1.000000000000000000e+00}" )),
                           std::forward_list<double>({-1.0, 0.0, 1.0}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<char> >( "{T, e, s, t,  , s, t, r, i, n, g}" )),
                           std::forward_list<char>({'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g'}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<std::string> >( "{Test, string}" )),
                           std::forward_list<std::string>({"Test","string"}) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<std::pair<int, int> > >( "{{0, 1}, {-1, 2}}" )),
                           (std::forward_list<std::pair<int,int> >({std::make_pair(0, 1), std::make_pair(-1, 2)})) );
  TEST_COMPARE_CONTAINERS( (Utility::fromString<std::forward_list<std::tuple<unsigned, double, long> > >( "{{0, 1.0, -100000}, {1, -1.00, 100001}}" )),
                           (std::forward_list<std::tuple<unsigned,double,long> >({std::make_tuple(0u, 1.0, -100000l), std::make_tuple(1u, -1.0, 100001l)})) );
}

//---------------------------------------------------------------------------//
// Check that a forward_list can be extracted from a stream
TEUCHOS_UNIT_TEST( ForwardList, fromStream )
{
  // Extract forward_list of short
  std::istringstream iss( "{-1, 2}" );

  {
    std::forward_list<short> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<short>({-1, 2}) );

    iss.str( "{-1, 2i, 2}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<short>({-1, 0, 1, 2}) );
  }

  // Extract list of unsigned short
  iss.str( "{0, 10, 100}" );
  iss.clear();

  {
    std::forward_list<unsigned short> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned short>({0, 10, 100}) );

    iss.str( "{0, 3i, 100}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned short>({0, 25, 50, 75, 100}) );
  }

  // Extract list of int
  iss.str( "{-11111, 0, 11111, 22222}" );
  iss.clear();

  {
    std::forward_list<int> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<int>({-11111, 0, 11111, 22222}) );

    iss.str( "{-11111, 2i, 22222}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<int>({-11111, 0, 11111, 22222}) );
  }

  // Extract list of unsigned int
  iss.str( "{0, 10, 100, 1000, 10000}" );
  iss.clear();

  {
    std::forward_list<unsigned> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned>({0, 10, 100, 1000, 10000}) );

    iss.str( "{0, 3i, 10000}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned>({0, 2500, 5000, 7500, 10000}) );
  }

  // Extract list of long
  iss.str( "{-11111, 0, 11111, 22222}" );
  iss.clear();

  {
    std::forward_list<long> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<long>({-11111, 0, 11111, 22222}) );

    iss.str( "{-11111, 2i, 22222}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<long>({-11111, 0, 11111, 22222}) );
  }

  // Extract list of unsigned long
  iss.str( "{0, 10, 100, 1000, 10000}" );
  iss.clear();

  {
    std::forward_list<unsigned long> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned long>({0, 10, 100, 1000, 10000}) );

    iss.str( "{0, 3i, 10000}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned long>({0, 2500, 5000, 7500, 10000}) );
  }

  // Extract list of long long
  iss.str( "{-1000000000, 0, 1000000000}" );
  iss.clear();

  {
    std::forward_list<long long> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<long long>({-1000000000, 0, 1000000000}) );

    iss.str( "{-1000000000, 1i, 1000000000}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<long long>({-1000000000, 0, 1000000000}) );
  }

  // Extract list of unsigned long long
  iss.str( "{0, 1000000000, 10000000000}" );
  iss.clear();

  {
    std::forward_list<unsigned long long> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned long long>({0, 1000000000, 10000000000}) );

    iss.str( "{0, 1i, 10000000000}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned long long>({0, 5000000000, 10000000000}) );
  }

  // Extract list of float
  iss.str( "{-1, 0.0, 1.000000000e+00}" );
  iss.clear();

  {
    std::forward_list<float> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<float>({-1.0f, 0.0f, 1.0f}) );

    iss.str( "{-1, 3i, 1.000000000e+00}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<float>({-1.0f, -0.5f, 0.0f, 0.5f, 1.0f}) );

    iss.str( "{1e-3, 2l, 1.0}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::forward_list<float>({1e-3f, 1e-2f, 1e-1f, 1.0f}),
                                      1e-7 );

    iss.str( "{1.0, 1l, 100, 3i, 200}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<float>({1.0f, 10.0f, 100.0f, 125.0f, 150.0f, 175.0f, 200.0f}) );

    iss.str( "{0.0, pi/2, 3Pi / 4, PI, 2*pi}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<float>({0.0f, (float)Utility::PhysicalConstants::pi/2, 3*(float)Utility::PhysicalConstants::pi/4, (float)Utility::PhysicalConstants::pi, 2*(float)Utility::PhysicalConstants::pi}) );

    iss.str( "{-pi, 3i, 0}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::forward_list<float>({-(float)Utility::PhysicalConstants::pi, -3*(float)Utility::PhysicalConstants::pi/4, -(float)Utility::PhysicalConstants::pi/2, -(float)Utility::PhysicalConstants::pi/4, 0.0}),
                                      1e-7 );

    iss.str( "{-inf, 0, Infinity}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<float>({-std::numeric_limits<float>::infinity(), 0.0f, std::numeric_limits<float>::infinity()}) );
  }

  // Extract list of double
  iss.str( "{-1, 0.0, 1.000000000000000000e+00}" );
  iss.clear();

  {
    std::forward_list<double> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<double>({-1.0, 0.0, 1.0}) );

    iss.str( "{-1, 3i, 1.000000000e+00}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<double>({-1.0, -0.5, 0.0, 0.5, 1.0}) );

    iss.str( "{1e-3, 2l, 1.0}" );
    iss.clear();

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::forward_list<double>({1e-3, 1e-2, 1e-1, 1.0}),
                                      1e-9 );

    iss.str( "{1.0, 1l, 100, 3i, 200}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::forward_list<double>({1.0, 10.0, 100.0, 125.0, 150.0, 175.0, 200.0}),
                                      1e-15 );

    iss.str( "{0.0, pi/2, 3Pi / 4, PI, 2*pi}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<double>({0.0, Utility::PhysicalConstants::pi/2, 3*Utility::PhysicalConstants::pi/4, Utility::PhysicalConstants::pi, 2*Utility::PhysicalConstants::pi}) );

    iss.str( "{-pi, 3i, 0}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<double>({-Utility::PhysicalConstants::pi, -3*Utility::PhysicalConstants::pi/4, -Utility::PhysicalConstants::pi/2, -Utility::PhysicalConstants::pi/4, 0.0}) );

    iss.str( "{-inf, 0, Infinity}" );
    iss.clear();

    Utility::fromStream( iss, test_list );
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<double>({-std::numeric_limits<float>::infinity(), 0.0, std::numeric_limits<float>::infinity()}) );
  }

  // Extract list of char
  iss.str( "{T, e, s, t,  , s, t, r, i, n, g}" );
  iss.clear();

  {
    std::forward_list<char> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<char>({'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g'}) );
  }

  // Extract list of string
  iss.str( "{Test, string}" );
  iss.clear();

  {
    std::list<std::string> test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<std::string>({"Test","string"}) );
  }

  // Extract list of pair
  iss.str( "{{0, 1}, {-1, 2}}" );
  iss.clear();

  {
    std::forward_list<std::pair<int,int> > test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, (std::forward_list<std::pair<int,int> >({std::make_pair(0, 1), std::make_pair(-1, 2)})) );
  }

  // Extract list of tuple
  iss.str( "{{0, 1.0, -100000}, {1, -1.00, 100001}}" );
  iss.clear();

  {
    std::forward_list<std::tuple<unsigned,double,long> > test_list;

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, (std::forward_list<std::tuple<unsigned,double,long> >({std::make_tuple(0u, 1.0, -100000l), std::make_tuple(1u, -1.0, 100001l)})) );
  }

  // Extract multiple lists of pairs
  iss.str( "{{0, 1}, {-1, 2}}, {{1, 0}, {2, -1}}" );
  iss.clear();

  {
    std::forward_list<std::pair<long,long> > test_list;

    Utility::fromStream( iss, test_list, "," );

    TEST_COMPARE_CONTAINERS( test_list, (std::forward_list<std::pair<long,long> >({std::make_pair(0l, 1l), std::make_pair(-1l, 2l)})) );

    Utility::moveInputStreamToNextElement( iss, ',', '}' );

    Utility::fromStream( iss, test_list );

    TEST_COMPARE_CONTAINERS( test_list, (std::forward_list<std::pair<long,long> >({std::make_pair(1l, 0l), std::make_pair(2l, -1l)})) );
  }

  // Extract list of lists
  iss.str( "{{-1, 3i, 1}, {1e-3, 2l, 1.0}, {-Inf, Inf}, {-pi/2, 1i, pi/2}}" );
  iss.clear();

  {
    std::forward_list<std::forward_list<double> > test_list;
    
    Utility::fromStream( iss, test_list );
    
    TEST_EQUALITY_CONST( std::distance( test_list.begin(), test_list.end() ), 4 );

    std::forward_list<std::forward_list<double> >::const_iterator
      test_list_it = test_list.begin();
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::forward_list<double>({-1.0, -0.5, 0.0, 0.5, 1.0})) );
    
    ++test_list_it;
    
    TEST_COMPARE_FLOATING_CONTAINERS( *test_list_it,
                                      (std::forward_list<double>({1e-3, 1e-2, 1e-1, 1.0})),
                                      1e-9 );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::forward_list<double>({-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()})) );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::forward_list<double>({-Utility::PhysicalConstants::pi/2, 0.0, Utility::PhysicalConstants::pi/2})) );
  }
}

//---------------------------------------------------------------------------//
// Check that a list can be placed in a stream
TEUCHOS_UNIT_TEST( List, ostream )
{
  std::ostringstream oss;

  oss << std::list<short>( {-1, 2} );
  
  TEST_EQUALITY_CONST( oss.str(), "{-1, 2}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<unsigned short>( {0, 10, 100} );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<int>( {-11111, 0, 11111, 22222} );
  
  TEST_EQUALITY_CONST( oss.str(), "{-11111, 0, 11111, 22222}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<unsigned int>( {0, 10, 100, 1000, 10000} );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100, 1000, 10000}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<long>( {-11111, 0, 11111, 22222} );
  
  TEST_EQUALITY_CONST( oss.str(), "{-11111, 0, 11111, 22222}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<unsigned long>( {0, 10, 100, 1000, 10000} );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100, 1000, 10000}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<long long>( {-10000000000, 0, 10000000000} );
  
  TEST_EQUALITY_CONST( oss.str(), "{-10000000000, 0, 10000000000}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<unsigned long long>( {0, 1000000000, 1000000000000} );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 1000000000, 1000000000000}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<float>( {-1.0f, 0.0f, 1.0f} );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{-1.000000000e+00, 0.000000000e+00, 1.000000000e+00}" );
  
  oss.str( "" );
  oss.clear();

  oss << std::list<double>( {-1.0, 0.0, 1.0} );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{-1.000000000000000000e+00, 0.000000000000000000e+00, 1.000000000000000000e+00}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<char>( {'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g' } );
  
  TEST_EQUALITY_CONST( oss.str(), "{T, e, s, t,  , s, t, r, i, n, g}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<std::string>( {"Test", "string"} );
  
  TEST_EQUALITY_CONST( oss.str(), "{Test, string}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<std::pair<int, int>>( {std::pair<int, int>({0, 1}), std::pair<int,int>({-1, 2})} );
  
  TEST_EQUALITY_CONST( oss.str(), "{{0, 1}, {-1, 2}}" );

  oss.str( "" );
  oss.clear();

  oss << std::list<std::tuple<unsigned, double, long>>( {std::make_tuple( 0u, 1.0, -100000 ), std::make_tuple( 1u, -1.0, 100001 )} );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{{0, 1.000000000000000000e+00, -100000}, {1, -1.000000000000000000e+00, 100001}}" );
}

//---------------------------------------------------------------------------//
// Check that a list can be extracted from a stream
TEUCHOS_UNIT_TEST( List, istream )
{
  // Extract list of short
  std::istringstream iss( "{-1, 2}" );

  {
    std::list<short> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<short>({-1, 2}) );

    iss.str( "{-1, 2i, 2}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<short>({-1, 0, 1, 2}) );
  }

  // Extract list of unsigned short
  iss.str( "{0, 10, 100}" );
  iss.clear();

  {
    std::list<unsigned short> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned short>({0, 10, 100}) );

    iss.str( "{0, 3i, 100}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned short>({0, 25, 50, 75, 100}) );
  }

  // Extract list of int
  iss.str( "{-11111, 0, 11111, 22222}" );
  iss.clear();

  {
    std::list<int> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<int>({-11111, 0, 11111, 22222}) );

    iss.str( "{-11111, 2i, 22222}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<int>({-11111, 0, 11111, 22222}) );
  }

  // Extract list of unsigned int
  iss.str( "{0, 10, 100, 1000, 10000}" );
  iss.clear();

  {
    std::list<unsigned> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned>({0, 10, 100, 1000, 10000}) );

    iss.str( "{0, 3i, 10000}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned>({0, 2500, 5000, 7500, 10000}) );
  }

  // Extract list of long
  iss.str( "{-11111, 0, 11111, 22222}" );
  iss.clear();

  {
    std::list<long> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<long>({-11111, 0, 11111, 22222}) );

    iss.str( "{-11111, 2i, 22222}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<long>({-11111, 0, 11111, 22222}) );
  }

  // Extract list of unsigned long
  iss.str( "{0, 10, 100, 1000, 10000}" );
  iss.clear();

  {
    std::list<unsigned long> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned long>({0, 10, 100, 1000, 10000}) );

    iss.str( "{0, 3i, 10000}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned long>({0, 2500, 5000, 7500, 10000}) );
  }

  // Extract list of long long
  iss.str( "{-1000000000, 0, 1000000000}" );
  iss.clear();

  {
    std::list<long long> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<long long>({-1000000000, 0, 1000000000}) );

    iss.str( "{-1000000000, 1i, 1000000000}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<long long>({-1000000000, 0, 1000000000}) );
  }

  // Extract list of unsigned long long
  iss.str( "{0, 1000000000, 10000000000}" );
  iss.clear();

  {
    std::list<unsigned long long> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned long long>({0, 1000000000, 10000000000}) );

    iss.str( "{0, 1i, 10000000000}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<unsigned long long>({0, 5000000000, 10000000000}) );
  }

  // Extract list of float
  iss.str( "{-1, 0.0, 1.000000000e+00}" );
  iss.clear();

  {
    std::list<float> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<float>({-1.0f, 0.0f, 1.0f}) );

    iss.str( "{-1, 3i, 1.000000000e+00}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<float>({-1.0f, -0.5f, 0.0f, 0.5f, 1.0f}) );

    iss.str( "{1e-3, 2l, 1.0}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::list<float>({1e-3f, 1e-2f, 1e-1f, 1.0f}),
                                      1e-7);

    iss.str( "{1.0, 1l, 100, 3i, 200}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<float>({1.0f, 10.0f, 100.0f, 125.0f, 150.0f, 175.0f, 200.0f}) );

    iss.str( "{0.0, pi/2, 3Pi / 4, PI, 2*pi}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<float>({0.0f, (float)Utility::PhysicalConstants::pi/2, 3*(float)Utility::PhysicalConstants::pi/4, (float)Utility::PhysicalConstants::pi, 2*(float)Utility::PhysicalConstants::pi}) );

    iss.str( "{-pi, 3i, 0}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::list<float>({-(float)Utility::PhysicalConstants::pi, -3*(float)Utility::PhysicalConstants::pi/4, -(float)Utility::PhysicalConstants::pi/2, -(float)Utility::PhysicalConstants::pi/4, 0.0}),
                                      1e-7 );

    iss.str( "{-inf, 0, Infinity}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<float>({-std::numeric_limits<float>::infinity(), 0.0f, std::numeric_limits<float>::infinity()}) );
  }

  // Extract list of double
  iss.str( "{-1, 0.0, 1.000000000000000000e+00}" );
  iss.clear();

  {
    std::list<double> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<double>({-1.0, 0.0, 1.0}) );

    iss.str( "{-1, 3i, 1.000000000e+00}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<double>({-1.0, -0.5, 0.0, 0.5, 1.0}) );

    iss.str( "{1e-3, 2l, 1.0}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::list<double>({1e-3, 1e-2, 1e-1, 1.0}),
                                      1e-9 );

    iss.str( "{1.0, 1l, 100, 3i, 200}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::list<double>({1.0, 10.0, 100.0, 125.0, 150.0, 175.0, 200.0}),
                                      1e-15 );

    iss.str( "{0.0, pi/2, 3Pi / 4, PI, 2*pi}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<double>({0.0, Utility::PhysicalConstants::pi/2, 3*Utility::PhysicalConstants::pi/4, Utility::PhysicalConstants::pi, 2*Utility::PhysicalConstants::pi}) );

    iss.str( "{-pi, 3i, 0}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<double>({-Utility::PhysicalConstants::pi, -3*Utility::PhysicalConstants::pi/4, -Utility::PhysicalConstants::pi/2, -Utility::PhysicalConstants::pi/4, 0.0}) );

    iss.str( "{-inf, 0, Infinity}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<double>({-std::numeric_limits<float>::infinity(), 0.0, std::numeric_limits<float>::infinity()}) );
  }

  // Extract list of char
  iss.str( "{T, e, s, t,  , s, t, r, i, n, g}" );
  iss.clear();

  {
    std::list<char> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<char>({'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g'}) );
  }

  // Extract list of string
  iss.str( "{Test, string}" );
  iss.clear();

  {
    std::list<std::string> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::list<std::string>({"Test","string"}) );
  }

  // Extract list of pair
  iss.str( "{{0, 1}, {-1, 2}}" );
  iss.clear();

  {
    std::list<std::pair<int,int> > test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, (std::list<std::pair<int,int> >({std::make_pair(0, 1), std::make_pair(-1, 2)})) );
  }

  // Extract list of tuple
  iss.str( "{{0, 1.0, -100000}, {1, -1.00, 100001}}" );
  iss.clear();

  {
    std::list<std::tuple<unsigned,double,long> > test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, (std::list<std::tuple<unsigned,double,long> >({std::make_tuple(0u, 1.0, -100000l), std::make_tuple(1u, -1.0, 100001l)})) );
  }

  // Extract multiple lists of pairs
  iss.str( "{{0, 1}, {-1, 2}}, {{1, 0}, {2, -1}}" );
  iss.clear();

  {
    std::list<std::pair<long,long> > test_list;

    Utility::fromStream( iss, test_list, "," );

    TEST_COMPARE_CONTAINERS( test_list, (std::list<std::pair<long,long> >({std::make_pair(0l, 1l), std::make_pair(-1l, 2l)})) );

    Utility::moveInputStreamToNextElement( iss, ',', '}' );

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, (std::list<std::pair<long,long> >({std::make_pair(1l, 0l), std::make_pair(2l, -1l)})) );
  }

  // Extract list of lists
  iss.str( "{{-1, 3i, 1}, {1e-3, 2l, 1.0}, {-Inf, Inf}, {-pi/2, 1i, pi/2}}" );
  iss.clear();

  {
    std::list<std::list<double> > test_list;
    
    iss >> test_list;
    
    TEST_EQUALITY_CONST( test_list.size(), 4 );

    std::list<std::list<double> >::const_iterator
      test_list_it = test_list.begin();
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-1.0, -0.5, 0.0, 0.5, 1.0})) );
    
    ++test_list_it;
    
    TEST_COMPARE_FLOATING_CONTAINERS( *test_list_it,
                                      (std::list<double>({1e-3, 1e-2, 1e-1, 1.0})),
                                      1e-9 );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()})) );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-Utility::PhysicalConstants::pi/2, 0.0, Utility::PhysicalConstants::pi/2})) );
  }
}

//---------------------------------------------------------------------------//
// Check that a forward_list can be placed in a stream
TEUCHOS_UNIT_TEST( ForwardList, ostream )
{
  std::ostringstream oss;

  oss << std::forward_list<short>( {-1, 2} );
  
  TEST_EQUALITY_CONST( oss.str(), "{-1, 2}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<unsigned short>( {0, 10, 100} );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<int>( {-11111, 0, 11111, 22222} );
  
  TEST_EQUALITY_CONST( oss.str(), "{-11111, 0, 11111, 22222}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<unsigned int>( {0, 10, 100, 1000, 10000} );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100, 1000, 10000}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<long>( {-11111, 0, 11111, 22222} );
  
  TEST_EQUALITY_CONST( oss.str(), "{-11111, 0, 11111, 22222}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<unsigned long>( {0, 10, 100, 1000, 10000} );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 10, 100, 1000, 10000}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<long long>( {-10000000000, 0, 10000000000} );
  
  TEST_EQUALITY_CONST( oss.str(), "{-10000000000, 0, 10000000000}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<unsigned long long>( {0, 1000000000, 1000000000000} );
  
  TEST_EQUALITY_CONST( oss.str(), "{0, 1000000000, 1000000000000}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<float>( {-1.0f, 0.0f, 1.0f} );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{-1.000000000e+00, 0.000000000e+00, 1.000000000e+00}" );
  
  oss.str( "" );
  oss.clear();

  oss << std::forward_list<double>( {-1.0, 0.0, 1.0} );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{-1.000000000000000000e+00, 0.000000000000000000e+00, 1.000000000000000000e+00}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<char>( {'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g' } );
  
  TEST_EQUALITY_CONST( oss.str(), "{T, e, s, t,  , s, t, r, i, n, g}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<std::string>( {"Test", "string"} );
  
  TEST_EQUALITY_CONST( oss.str(), "{Test, string}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<std::pair<int, int>>( {std::pair<int, int>({0, 1}), std::pair<int,int>({-1, 2})} );
  
  TEST_EQUALITY_CONST( oss.str(), "{{0, 1}, {-1, 2}}" );

  oss.str( "" );
  oss.clear();

  oss << std::forward_list<std::tuple<unsigned, double, long>>( {std::make_tuple( 0u, 1.0, -100000 ), std::make_tuple( 1u, -1.0, 100001 )} );
  
  TEST_EQUALITY_CONST( oss.str(),
                       "{{0, 1.000000000000000000e+00, -100000}, {1, -1.000000000000000000e+00, 100001}}" );
}

//---------------------------------------------------------------------------//
// Check that a forward_list can be extracted from a stream
TEUCHOS_UNIT_TEST( ForwardList, istream )
{
  // Extract forward_list of short
  std::istringstream iss( "{-1, 2}" );

  {
    std::forward_list<short> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<short>({-1, 2}) );

    iss.str( "{-1, 2i, 2}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<short>({-1, 0, 1, 2}) );
  }

  // Extract list of unsigned short
  iss.str( "{0, 10, 100}" );
  iss.clear();

  {
    std::forward_list<unsigned short> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned short>({0, 10, 100}) );

    iss.str( "{0, 3i, 100}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned short>({0, 25, 50, 75, 100}) );
  }

  // Extract list of int
  iss.str( "{-11111, 0, 11111, 22222}" );
  iss.clear();

  {
    std::forward_list<int> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<int>({-11111, 0, 11111, 22222}) );

    iss.str( "{-11111, 2i, 22222}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<int>({-11111, 0, 11111, 22222}) );
  }

  // Extract list of unsigned int
  iss.str( "{0, 10, 100, 1000, 10000}" );
  iss.clear();

  {
    std::forward_list<unsigned> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned>({0, 10, 100, 1000, 10000}) );

    iss.str( "{0, 3i, 10000}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned>({0, 2500, 5000, 7500, 10000}) );
  }

  // Extract list of long
  iss.str( "{-11111, 0, 11111, 22222}" );
  iss.clear();

  {
    std::forward_list<long> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<long>({-11111, 0, 11111, 22222}) );

    iss.str( "{-11111, 2i, 22222}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<long>({-11111, 0, 11111, 22222}) );
  }

  // Extract list of unsigned long
  iss.str( "{0, 10, 100, 1000, 10000}" );
  iss.clear();

  {
    std::forward_list<unsigned long> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned long>({0, 10, 100, 1000, 10000}) );

    iss.str( "{0, 3i, 10000}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned long>({0, 2500, 5000, 7500, 10000}) );
  }

  // Extract list of long long
  iss.str( "{-1000000000, 0, 1000000000}" );
  iss.clear();

  {
    std::forward_list<long long> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<long long>({-1000000000, 0, 1000000000}) );

    iss.str( "{-1000000000, 1i, 1000000000}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<long long>({-1000000000, 0, 1000000000}) );
  }

  // Extract list of unsigned long long
  iss.str( "{0, 1000000000, 10000000000}" );
  iss.clear();

  {
    std::forward_list<unsigned long long> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned long long>({0, 1000000000, 10000000000}) );

    iss.str( "{0, 1i, 10000000000}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<unsigned long long>({0, 5000000000, 10000000000}) );
  }

  // Extract list of float
  iss.str( "{-1, 0.0, 1.000000000e+00}" );
  iss.clear();

  {
    std::forward_list<float> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<float>({-1.0f, 0.0f, 1.0f}) );

    iss.str( "{-1, 3i, 1.000000000e+00}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<float>({-1.0f, -0.5f, 0.0f, 0.5f, 1.0f}) );

    iss.str( "{1e-3, 2l, 1.0}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::forward_list<float>({1e-3f, 1e-2f, 1e-1f, 1.0f}),
                                      1e-7 );

    iss.str( "{1.0, 1l, 100, 3i, 200}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<float>({1.0f, 10.0f, 100.0f, 125.0f, 150.0f, 175.0f, 200.0f}) );

    iss.str( "{0.0, pi/2, 3Pi / 4, PI, 2*pi}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<float>({0.0f, (float)Utility::PhysicalConstants::pi/2, 3*(float)Utility::PhysicalConstants::pi/4, (float)Utility::PhysicalConstants::pi, 2*(float)Utility::PhysicalConstants::pi}) );

    iss.str( "{-pi, 3i, 0}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::forward_list<float>({-(float)Utility::PhysicalConstants::pi, -3*(float)Utility::PhysicalConstants::pi/4, -(float)Utility::PhysicalConstants::pi/2, -(float)Utility::PhysicalConstants::pi/4, 0.0}),
                                      1e-7 );

    iss.str( "{-inf, 0, Infinity}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<float>({-std::numeric_limits<float>::infinity(), 0.0f, std::numeric_limits<float>::infinity()}) );
  }

  // Extract list of double
  iss.str( "{-1, 0.0, 1.000000000000000000e+00}" );
  iss.clear();

  {
    std::forward_list<double> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<double>({-1.0, 0.0, 1.0}) );

    iss.str( "{-1, 3i, 1.000000000e+00}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<double>({-1.0, -0.5, 0.0, 0.5, 1.0}) );

    iss.str( "{1e-3, 2l, 1.0}" );
    iss.clear();

    iss >> test_list;

    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::forward_list<double>({1e-3, 1e-2, 1e-1, 1.0}),
                                      1e-9 );

    iss.str( "{1.0, 1l, 100, 3i, 200}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_FLOATING_CONTAINERS( test_list,
                                      std::forward_list<double>({1.0, 10.0, 100.0, 125.0, 150.0, 175.0, 200.0}),
                                      1e-15 );

    iss.str( "{0.0, pi/2, 3Pi / 4, PI, 2*pi}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<double>({0.0, Utility::PhysicalConstants::pi/2, 3*Utility::PhysicalConstants::pi/4, Utility::PhysicalConstants::pi, 2*Utility::PhysicalConstants::pi}) );

    iss.str( "{-pi, 3i, 0}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<double>({-Utility::PhysicalConstants::pi, -3*Utility::PhysicalConstants::pi/4, -Utility::PhysicalConstants::pi/2, -Utility::PhysicalConstants::pi/4, 0.0}) );

    iss.str( "{-inf, 0, Infinity}" );
    iss.clear();

    iss >> test_list;
    
    TEST_COMPARE_CONTAINERS( test_list,
                             std::forward_list<double>({-std::numeric_limits<float>::infinity(), 0.0, std::numeric_limits<float>::infinity()}) );
  }

  // Extract list of char
  iss.str( "{T, e, s, t,  , s, t, r, i, n, g}" );
  iss.clear();

  {
    std::forward_list<char> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<char>({'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g'}) );
  }

  // Extract list of string
  iss.str( "{Test, string}" );
  iss.clear();

  {
    std::list<std::string> test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, std::forward_list<std::string>({"Test","string"}) );
  }

  // Extract list of pair
  iss.str( "{{0, 1}, {-1, 2}}" );
  iss.clear();

  {
    std::forward_list<std::pair<int,int> > test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, (std::forward_list<std::pair<int,int> >({std::make_pair(0, 1), std::make_pair(-1, 2)})) );
  }

  // Extract list of tuple
  iss.str( "{{0, 1.0, -100000}, {1, -1.00, 100001}}" );
  iss.clear();

  {
    std::forward_list<std::tuple<unsigned,double,long> > test_list;

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, (std::forward_list<std::tuple<unsigned,double,long> >({std::make_tuple(0u, 1.0, -100000l), std::make_tuple(1u, -1.0, 100001l)})) );
  }

  // Extract multiple lists of pairs
  iss.str( "{{0, 1}, {-1, 2}}, {{1, 0}, {2, -1}}" );
  iss.clear();

  {
    std::forward_list<std::pair<long,long> > test_list;

    Utility::fromStream( iss, test_list, "," );

    TEST_COMPARE_CONTAINERS( test_list, (std::forward_list<std::pair<long,long> >({std::make_pair(0l, 1l), std::make_pair(-1l, 2l)})) );

    Utility::moveInputStreamToNextElement( iss, ',', '}' );

    iss >> test_list;

    TEST_COMPARE_CONTAINERS( test_list, (std::forward_list<std::pair<long,long> >({std::make_pair(1l, 0l), std::make_pair(2l, -1l)})) );
  }

  // Extract list of lists
  iss.str( "{{-1, 3i, 1}, {1e-3, 2l, 1.0}, {-Inf, Inf}, {-pi/2, 1i, pi/2}}" );
  iss.clear();

  {
    std::forward_list<std::forward_list<double> > test_list;
    
    iss >> test_list;
    
    TEST_EQUALITY_CONST( std::distance( test_list.begin(), test_list.end() ), 4 );

    std::forward_list<std::forward_list<double> >::const_iterator
      test_list_it = test_list.begin();
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::forward_list<double>({-1.0, -0.5, 0.0, 0.5, 1.0})) );
    
    ++test_list_it;
    
    TEST_COMPARE_FLOATING_CONTAINERS( *test_list_it,
                                      (std::forward_list<double>({1e-3, 1e-2, 1e-1, 1.0})),
                                      1e-9 );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::forward_list<double>({-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()})) );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::forward_list<double>({-Utility::PhysicalConstants::pi/2, 0.0, Utility::PhysicalConstants::pi/2})) );
  }
}

//---------------------------------------------------------------------------//
// Check that a list can mimic stream behavior
TEUCHOS_UNIT_TEST( List, stream_mimic )
{
  {
    std::list<short> test_list;

    test_list << -1 << 0 << 1;

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<short>( {-1, 0, 1} ) );

    short element0, element1, element2;

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, -1 );
    TEST_EQUALITY_CONST( element1, 0 );
    TEST_EQUALITY_CONST( element2, 1 );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "-1" << 0 << "1";

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<short>( {-1, 0, 1} ) );

    std::string element0_str, element1_str, element2_str;

    test_list >> element2_str >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "-1" );
    TEST_EQUALITY_CONST( element1_str, "0" );
    TEST_EQUALITY_CONST( element2_str, "1" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }

  {
    std::list<unsigned short> test_list;

    test_list << 0 << 10 << 100;

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<unsigned short>( {0, 10, 100} ) );

    unsigned short element0, element1, element2;

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, 0 );
    TEST_EQUALITY_CONST( element1, 10 );
    TEST_EQUALITY_CONST( element2, 100 );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "0" << 10 << "100";

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<unsigned short>( {0, 10, 100} ) );

    std::string element0_str, element1_str, element2_str;

    test_list >> element2_str >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "0" );
    TEST_EQUALITY_CONST( element1_str, "10" );
    TEST_EQUALITY_CONST( element2_str, "100" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }

  {
    std::list<int> test_list;

    test_list << -11111 << 0 << 11111;

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<int>( {-11111, 0, 11111} ) );

    int element0, element1, element2;

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, -11111 );
    TEST_EQUALITY_CONST( element1, 0 );
    TEST_EQUALITY_CONST( element2, 11111 );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "-11111" << 0 << "11111";

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<int>( {-11111, 0, 11111} ) );

    std::string element0_str, element1_str, element2_str;

    test_list >> element2_str >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "-11111" );
    TEST_EQUALITY_CONST( element1_str, "0" );
    TEST_EQUALITY_CONST( element2_str, "11111" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }

  {
    std::list<unsigned int> test_list;

    test_list << 0 << 1000 << 1000000;

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<unsigned int>( {0, 1000, 1000000} ) );

    unsigned int element0, element1, element2;

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, 0 );
    TEST_EQUALITY_CONST( element1, 1000 );
    TEST_EQUALITY_CONST( element2, 1000000 );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "0" << 1000 << "1000000";

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<unsigned int>( {0, 1000, 1000000} ) );

    std::string element0_str, element1_str, element2_str;

    test_list >> element2_str >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "0" );
    TEST_EQUALITY_CONST( element1_str, "1000" );
    TEST_EQUALITY_CONST( element2_str, "1000000" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }

  {
    std::list<long> test_list;

    test_list << -11111l << 0l << 11111l;

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<long>( {-11111l, 0l, 11111l} ) );

    long element0, element1, element2;

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, -11111l );
    TEST_EQUALITY_CONST( element1, 0l );
    TEST_EQUALITY_CONST( element2, 11111l );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "-11111" << 0 << "11111";

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<long>( {-11111l, 0l, 11111l} ) );

    std::string element0_str, element1_str, element2_str;

    test_list >> element2_str >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "-11111" );
    TEST_EQUALITY_CONST( element1_str, "0" );
    TEST_EQUALITY_CONST( element2_str, "11111" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }

  {
    std::list<unsigned long> test_list;

    test_list << 0ul << 1000ul << 1000000ul;

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<unsigned long>( {0ul, 1000ul, 1000000ul} ) );

    unsigned long element0, element1, element2;

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, 0ul );
    TEST_EQUALITY_CONST( element1, 1000ul );
    TEST_EQUALITY_CONST( element2, 1000000ul );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "0" << 1000 << "1000000";

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<unsigned long>( {0ul, 1000ul, 1000000ul} ) );

    std::string element0_str, element1_str, element2_str;

    test_list >> element2_str >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "0" );
    TEST_EQUALITY_CONST( element1_str, "1000" );
    TEST_EQUALITY_CONST( element2_str, "1000000" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }
  
  {
    std::list<long long> test_list;

    test_list << -10000000000ll << 0ll << 10000000000ll;

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<long long>( {-10000000000ll, 0ll, 10000000000ll} ) );

    long long element0, element1, element2;

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, -10000000000ll );
    TEST_EQUALITY_CONST( element1, 0ll );
    TEST_EQUALITY_CONST( element2, 10000000000ll );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "-10000000000" << 0 << "10000000000";

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<long long>( {-10000000000ll, 0ll, 10000000000ll} ) );

    std::string element0_str, element1_str, element2_str;

    test_list >> element2_str >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "-10000000000" );
    TEST_EQUALITY_CONST( element1_str, "0" );
    TEST_EQUALITY_CONST( element2_str, "10000000000" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }

  {
    std::list<unsigned long long> test_list;

    test_list << 0ull << 10000000000ull << 100000000000ull;

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<unsigned long long>( {0ull, 10000000000ull, 100000000000ull} ) );

    unsigned long long element0, element1, element2;

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, 0ull );
    TEST_EQUALITY_CONST( element1, 10000000000ull );
    TEST_EQUALITY_CONST( element2, 100000000000ull );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "0" << 10000000000ull << "100000000000";

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<unsigned long long>( {0ull, 10000000000ull, 100000000000ull} ) );

    std::string element0_str, element1_str, element2_str;

    test_list >> element2_str >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "0" );
    TEST_EQUALITY_CONST( element1_str, "10000000000" );
    TEST_EQUALITY_CONST( element2_str, "100000000000" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }

  {
    std::list<float> test_list;

    test_list << -1.0f << 0.0f << 1.0f;

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<float>( {-1.0f, 0.0f, 1.0f} ) );

    float element0, element1, element2;

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, -1.0f );
    TEST_EQUALITY_CONST( element1, 0.0f );
    TEST_EQUALITY_CONST( element2, 1.0f );

    test_list << "-2pi/3" << "0" << "Inf";

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, -2*(float)Utility::PhysicalConstants::pi/3 );
    TEST_EQUALITY_CONST( element1, 0.0f );
    TEST_EQUALITY_CONST( element2, std::numeric_limits<float>::infinity() );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << -1.0f << 0.0f << 1.0f;

    std::string element0_str, element1_str, element2_str;

    test_list >> element2_str >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "-1.000000000e+00" );
    TEST_EQUALITY_CONST( element1_str, "0.000000000e+00" );
    TEST_EQUALITY_CONST( element2_str, "1.000000000e+00" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }

  {
    std::list<double> test_list;

    test_list << -1.0 << 0.0 << 1.0;

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<double>( {-1.0, 0.0, 1.0} ) );

    double element0, element1, element2;

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, -1.0 );
    TEST_EQUALITY_CONST( element1, 0.0 );
    TEST_EQUALITY_CONST( element2, 1.0 );

    test_list << "-infty" << "0" << "pi/3";

    test_list >> element2 >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, -std::numeric_limits<double>::infinity() );
    TEST_EQUALITY_CONST( element1, 0.0 );
    TEST_EQUALITY_CONST( element2, Utility::PhysicalConstants::pi/3 );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << -1.0 << 0.0 << 1.0;

    std::string element0_str, element1_str, element2_str;

    test_list >> element2_str >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "-1.000000000000000000e+00" );
    TEST_EQUALITY_CONST( element1_str, "0.000000000000000000e+00" );
    TEST_EQUALITY_CONST( element2_str, "1.000000000000000000e+00" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }

  {
    std::list<char> test_list;

    test_list << 'T' << 'e' << 's' << 't' << ' '
                << 's' << 't' << 'r' << 'i' << 'n' << 'g';

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<char>( {'T', 'e', 's', 't', ' ', 's', 't', 'r', 'i', 'n', 'g'} ) );

    std::string elements( "           " );

    test_list >> elements[10] >> elements[9] >> elements[8] >> elements[7]
                >> elements[6] >> elements[5] >> elements[4] >> elements[3]
                >> elements[2] >> elements[1] >> elements[0];

    TEST_EQUALITY_CONST( elements, "Test string" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }
  
  {
    std::list<std::string> test_list;

    test_list << "Test" << "string";

    TEST_COMPARE_CONTAINERS( test_list,
                             std::list<std::string>( {"Test", "string"} ) );

    std::string element0, element1;

    test_list >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, "Test" );
    TEST_EQUALITY_CONST( element1, "string" );
  }

  {
    std::list<std::pair<int, int> > test_list;

    test_list << std::make_pair( 0, 1 ) << std::make_pair( -1, 2 );

    TEST_COMPARE_CONTAINERS( test_list,
                             (std::list<std::pair<int, int> >( {std::pair<int, int>({0, 1}), std::pair<int,int>({-1, 2})} )) );

    std::pair<int,int> element0, element1;

    test_list >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, std::make_pair( 0, 1 ) );
    TEST_EQUALITY_CONST( element1, std::make_pair( -1, 2 ) );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "{0, 1}" << "{-1, 2}";

    test_list >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, std::make_pair( 0, 1 ) );
    TEST_EQUALITY_CONST( element1, std::make_pair( -1, 2 ) );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "{0, 1}" << std::make_pair( -1, 2 );

    std::string element0_str, element1_str;

    test_list >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "{0, 1}" );
    TEST_EQUALITY_CONST( element1_str, "{-1, 2}" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }

  {
    std::list<std::tuple<unsigned, double, long> > test_list;

    test_list << std::make_tuple( 0u, 1.0, -100000l )
                << std::make_tuple( 1u, -1.0, 100001l );

    TEST_COMPARE_CONTAINERS( test_list,
                             (std::list<std::tuple<unsigned, double, long> >( {std::make_tuple( 0u, 1.0, -100000 ), std::make_tuple( 1u, -1.0, 100001 )} )) );

    std::tuple<unsigned, double, long> element0, element1;

    test_list >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, std::make_tuple( 0u, 1.0, -100000l ) );
    TEST_EQUALITY_CONST( element1, std::make_tuple( 1u, -1.0, 100001l ) );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "{0, 1, -100000}" << "{1, -1, 100001}";

    test_list >> element1 >> element0;

    TEST_EQUALITY_CONST( element0, std::make_tuple( 0u, 1.0, -100000l ) );
    TEST_EQUALITY_CONST( element1, std::make_tuple( 1u, -1.0, 100001l ) );
    TEST_EQUALITY_CONST( test_list.size(), 0 );

    test_list << "{0, 1, -100000}" << std::make_tuple( 1u, -1.0, 100001l );

    std::string element0_str, element1_str;

    test_list >> element1_str >> element0_str;

    TEST_EQUALITY_CONST( element0_str, "{0, 1.000000000000000000e+00, -100000}" );
    TEST_EQUALITY_CONST( element1_str, "{1, -1.000000000000000000e+00, 100001}" );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }

  {
    std::list<std::list<double> > test_list;

    test_list << std::list<double>( {-1, -0.5, 0.0, 0.5, 1.0} )
              << std::list<double>( {1e-3, 1e-2, 1e-1, 1.0} )
              << std::list<double>( {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()} )
              << std::list<double>( {-Utility::PhysicalConstants::pi/2, 0.0, Utility::PhysicalConstants::pi/2} );

    TEST_EQUALITY_CONST( test_list.size(), 4 );

    std::list<std::list<double> >::const_iterator
      test_list_it = test_list.begin();
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-1.0, -0.5, 0.0, 0.5, 1.0})) );

    ++test_list_it;
    
    TEST_COMPARE_FLOATING_CONTAINERS( *test_list_it,
                                      (std::list<double>({1e-3, 1e-2, 1e-1, 1.0})),
                                      1e-9 );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()})) );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-Utility::PhysicalConstants::pi/2, 0.0, Utility::PhysicalConstants::pi/2})) );

    std::list<double> element0, element1, element2, element3;
    
    test_list >> element3 >> element2 >> element1 >> element0;
    
    TEST_COMPARE_CONTAINERS( element0,
                             (std::list<double>({-1.0, -0.5, 0.0, 0.5, 1.0})) );
    TEST_COMPARE_FLOATING_CONTAINERS( element1,
                                      (std::list<double>({1e-3, 1e-2, 1e-1, 1.0})),
                                      1e-9 );
    TEST_COMPARE_CONTAINERS( element2,
                             (std::list<double>({-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()})) );
    TEST_COMPARE_CONTAINERS( element3,
                             (std::list<double>({-Utility::PhysicalConstants::pi/2, 0.0, Utility::PhysicalConstants::pi/2})) );

    test_list << "{-1, 3i, 1}" << "{1e-3, 2l, 1.0}"
                << "{-Inf, Inf}" << "{-pi/2, 1i, pi/2}";

    TEST_EQUALITY_CONST( test_list.size(), 4 );

    test_list_it = test_list.begin();
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-1.0, -0.5, 0.0, 0.5, 1.0})) );

    ++test_list_it;
    
    TEST_COMPARE_FLOATING_CONTAINERS( *test_list_it,
                                      (std::list<double>({1e-3, 1e-2, 1e-1, 1.0})),
                                      1e-9 );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()})) );

    ++test_list_it;
    
    TEST_COMPARE_CONTAINERS( *test_list_it,
                             (std::list<double>({-Utility::PhysicalConstants::pi/2, 0.0, Utility::PhysicalConstants::pi/2})) );

    std::string element0_str, element1_str, element2_str, element3_str;

    test_list >> element3_str >> element2_str
                >> element1_str >> element0_str;

    TEST_ASSERT( element0_str.size() > 0 );
    TEST_ASSERT( element1_str.size() > 0 );
    TEST_ASSERT( element2_str.size() > 0 );
    TEST_ASSERT( element3_str.size() > 0 );
    TEST_EQUALITY_CONST( test_list.size(), 0 );
  }
}

//---------------------------------------------------------------------------//
// end tstList.cpp
//---------------------------------------------------------------------------//
