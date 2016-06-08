#ifndef POISSON_PROBLEMS_HH
#define POISSON_PROBLEMS_HH

#include <cassert>
#include <cmath>

#include "temporalprobleminterface.hh"

template <class FunctionSpace>
class DiffusionProblem : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;

  DiffusionProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider )
  {}

  double pow( double x, int a ) const
  {
    if( a > 0 )
      return pow( x, a-1 ) * x;
    if( a < 0 )
      return 1.0 / pow( x, a );
    return 1.0;
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& ret) const
  {
    ret = RangeType(0);

#if TEST_CASE && !FIXED_MESH
    // define evolution of surface
    const double at = 1.0 + 0.25 * sin( time() );
    const double apt = 0.25 * cos( time() );

    // calculated surface parameters
    const double divGammaV = 0.5 * at * apt * ( x[1]*x[1] + x[2]*x[2] ) / ( x[0]*x[0] + at*at * ( x[1]*x[1] + x[2]*x[2] ) );
    const double N1 = 1/at * x[0] / sqrt( x[0]*x[0] / (at*at) + x[1]*x[1] + x[2]*x[2] );
    const double N2 = x[1] / sqrt( x[0]*x[0] / (at*at) + x[1]*x[1] + x[2]*x[2] );
    const double H = ( 2.0 * x[0] * x[0] + at * ( 1 + at ) * ( x[1]*x[1] + x[2]*x[2] ) )
      / ( sqrt( x[0]*x[0] / (at*at) + x[1]*x[1] + x[2]*x[2] ) * ( x[0]*x[0] + at*at * ( x[1]*x[1] + x[2]*x[2] ) ) );


    // calculate solution and derivatives
    const double ux = sin( time() ) * x[0] * x[1];
    const double mdux = ( cos( time() ) + 0.5 * sin( time() ) * apt / at ) * x[0] * x[1];
    const double mlapux = sin( time() ) * (  2.0 * N1 * N2 + H * ( x[1] * N1 + x[0] * N2 ) );

    // construct solution
    ret = mdux + divGammaV * ux + mlapux;
#endif
  }

#if TEST_CASE
  virtual void D( RangeType &DD ) const
  {
    DD = RangeType(1);
  }
#endif

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
#if TEST_CASE
#if FIXED_MESH
    phi = exp( - 6.0 * time() ) * x[0] * x[1];
#else
    phi = sin( time() ) * x[0] * x[1];
#endif
#else
    // no exact solution - method needed for initial data
    phi = RangeType(1);
#endif
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
#if TEST_CASE
#if FIXED_MESH
    // no exact solution
    ret[ 0 ][ 0 ] = exp( -6.0 * time() ) * x[1] * ( 1.0 - 2.0 * x[0] * x[0] );
    ret[ 0 ][ 1 ] = exp( -6.0 * time() ) * x[0] * ( 1.0 - 2.0 * x[1] * x[1] );
    ret[ 0 ][ 2 ] = exp( -6.0 * time() ) * -2.0 * x[0] * x[1] * x[2];
#else
    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = sin( time() ) * x[1];
    grad[ 0 ][ 1 ] = sin( time() ) * x[0];
    grad[ 0 ][ 2 ] = 0.0;

    const double at = 1.0 + 0.25 * sin( time() );
    DomainType nu;
    nu[ 0 ] = 2.0 * x[ 0 ] / at;
    nu[ 1 ] = 2.0 * x[ 1 ];
    nu[ 2 ] = 2.0 * x[ 2 ];
    nu /= nu.two_norm();

    double dot = 0;
    for( unsigned int i = 0; i < nu.size(); ++i )
      {
	dot += nu[ i ] * grad[ 0 ][ i ];
      }

    for( unsigned int i = 0; i < nu.size(); ++i )
      {
	ret[ 0 ][ i ] = grad[ 0 ][ i ] - dot * nu[ i ];
      }
#endif
#else
    // no exact solution
    ret = JacobianRangeType(0);
#endif
  }
};

template <class FunctionSpace>
class FrapProblem : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;

  FrapProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider ),
      c0( Dune::Fem::Parameter::getValue< double >("diffusion.c0", 1 ) ),
      cx( Dune::Fem::Parameter::getValue< double >("diffusion.frap.cx", 0.25 ) ),
      cy( Dune::Fem::Parameter::getValue< double >("diffusion.frap.cy", 0.25 ) ),
      cz( Dune::Fem::Parameter::getValue< double >("diffusion.frap.cz", 0.25 ) ),
      r( Dune::Fem::Parameter::getValue< double >("diffusion.frap.r", 0.25 ) )
  {}

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& ret) const
  {
    ret = RangeType(0);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
    phi = RangeType( c0 );

    if( region( x ) )
      phi = RangeType(0);
  }

  //! region method for helping with recovery tests
  virtual bool region( const DomainType& x ) const
  {
    DomainType c = { cx, cy, cz };
    if( ( c - x ).two_norm() < r )
      return true;

    return false;
  }

private:
  const double c0;
  const double cx, cy, cz;
  const double r;
};

template <class FunctionSpace>
class LigandProblem : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;

  LigandProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider ),
      c0_( Dune::Fem::Parameter::getValue< double >("diffusion.c0", 1 ) ),
      k_( Dune::Fem::Parameter::getValue< double >("diffusion.ligand.k", 1 ) ),
      d_( Dune::Fem::Parameter::getValue< int >("diffusion.ligand.dir", 1 ) ),
      cx_( Dune::Fem::Parameter::getValue< double >("diffusion.ligand.cx", 0.25 ) ),
      cy_( Dune::Fem::Parameter::getValue< double >("diffusion.ligand.cy", 0.25 ) ),
      cz_( Dune::Fem::Parameter::getValue< double >("diffusion.ligand.cz", 0.25 ) ),
      r_( Dune::Fem::Parameter::getValue< double >("diffusion.ligand.r", 0.25 ) )
  {}

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& ret) const
  {
    ret = RangeType(0);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
    phi = RangeType( c0_ );

    if( region( x ) )
      phi = RangeType( 0 );
  }

  //! mass coefficient (default = 0)
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = k_ * x[ d_ ];
  }

  //! region method for helping with recovery tests
  virtual bool region( const DomainType& x ) const
  {
    DomainType c = { cx_, cy_, cz_ };
    if( ( c - x ).two_norm() > r_ )
      return true;

    return false;
  }

private:
  const double c0_;
  const double k_;
  const int d_;
  const double cx_, cy_, cz_;
  const double r_;
};

#endif // #ifndef POISSON_HH
