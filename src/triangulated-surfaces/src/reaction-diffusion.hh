#ifndef REACTIONDIFFUSION_HH
#define REACTIONDIFFUSION_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>

#include "temporalprobleminterface.hh"
#include "model.hh"

template < class FunctionSpace >
class KochMeinhardtProblem : public TemporalProblemInterface < FunctionSpace >
{
 typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;
  using BaseType :: deltaT ;

  KochMeinhardtProblem( const TimeProviderType &timeProvider )
    : BaseType( timeProvider ),
      // initial values
      a0 ( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.a0" ) ),
      s0 ( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.s0" ) ),
      var0( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.initialvariation" ) ),
      // diffusion coefficients
      D1 ( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.D1" ) ),
      D2 ( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.D2" ) ),
      // reaction coefficients
      rho1 ( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.rho1" ) ),
      rho2  ( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.rho2" ) ),
      k ( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.k" ) ),
      mu1 ( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.mu1" ) ),
      mu2 ( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.mu2" ) ),
      sigma1 ( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.sigma1" ) ),
      sigma2 ( Dune::Fem::Parameter::getValue< double >( "KochMeinhardt.sigma2" ) )
  {}

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
    phi = RangeType(0);

#if TEST_CASE
    phi[ 0 ] = exp( -6.0 * D1 * time() ) * x[0] * x[1];
    phi[ 1 ] = exp( -6.0 * D2 * time() ) * x[0] * x[1];
#else
    phi[ 0 ] = a0 * ( 1.0 + var0 * ( ( double) rand() / RAND_MAX - 0.5 ) );
    phi[ 1 ] = s0 * ( 1.0 + var0 * ( ( double) rand() / RAND_MAX - 0.5 ) );
    // phi[ 0 ] = a0 * ( 1.0 + var0 * sin( 10.0 * x[0] ) );
    // phi[ 1 ] = s0 * ( 1.0 + var0 * sin( 10.0 * x[1] ) );
#endif
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& value) const
  {
#if TEST_CASE
    const double a = exp( -6.0 * D1 * time() ) * x[0] * x[1];
    const double s = exp( -6.0 * D2 * time() ) * x[0] * x[1];

    value[ 0 ] = - rho1 * a*a*s / ( 1 + k * a*a ) + mu1 * a;
    value[ 1 ] = + rho2 * a*a*s / ( 1 + k * a*a )- sigma2;
#else
    value = 0;
#endif
  }

  //! diffusion coefficient (default = Id)
  virtual void D( RangeType& D ) const
  {
    D = 0;

    // c1 diffusion coefficient
    D[ 0 ] = D1;
    // c2 diffusion coefficient
    D[ 1 ] = D2;
  }

  //! non linear term in reaction diffusion system
  virtual void nonLinearity( const RangeType &uBar,
			     const RangeType &u,
			     RangeType &ret,
			     const bool implicit ) const
  {
    /**
     *  full non-linearity is
     *   f_1(a,s) = rho1 * a*a*s / ( 1 + k * a*a ) - mu1 * a
     *   f_2(a,s) = -rho2 * a*a*s / ( 1 + k * a*a )+ sigma2
     *
     *  this is linearized as
     *   f(a,s) approx f(a_0, s_0) + f_a(a_0,s_0) ( a - a_0 )
     *                  + f_s(a_0,s_0) ( s - s_0 )
     *  that is
     *   f_1(a,s) approx
     *     rho1 * ( aold*aold*s / ( 1 + aold*aold*k )
     *              + 2 * a * aold*sold / ( 1 + aold*aold*k )^2
     *              - 2 * aold*aold*sold / ( 1 + aold*aold*k )^2 )
     *        -mu1 * a
     *   f_2(a,s) approx
     *     -rho2 * ( aold*aold*s / ( 1 + aold*aold*k )
     *              + 2 * a * aold*sold / ( 1 + aold*aold*k )^2
     *              - 2 * aold*aold*sold / ( 1 + aold*aold*k )^2 )
     *        +sigma2
     *
     *  here this is implemented as
     *  implicit:
     *    ret[0] = -rho1 * ( u[1]*uBar[0]*uBar[0] / ( 1 + uBar[0]*uBar[0]*k )
     *                      + 2.0 * u[0]*uBar[0]*uBar[1] / ( (1 + uBar[0]*uBar[0]*k)*(1 + uBar[0]*uBar[0]*k) ) )
     *              - mu1 * u[0];
     *    ret[1] = rho2 * ( u[1]*uBar[0]*uBar[0] / ( 1 + uBar[0]*uBar[0]*k )
     *                      + 2.0 * u[0]*uBar[0]*uBar[1] / ( (1 + uBar[0]*uBar[0]*k)*(1 + uBar[0]*uBar[0]*k) ) );
     *  explicit:
     *    ret[0] = rho1 * ( -2.0 * uBar[0]*uBar[0]*uBar[1] / ( (1 + uBar[0]*uBar[0]*k)*(1 + uBar[0]*uBar[0]*k) ) );
     *    ret[1] = -rho2 * ( -2.0 * uBar[0]*uBar[0]*uBar[1] / ( (1 + uBar[0]*uBar[0]*k)*(1 + uBar[0]*uBar[0]*k) ) )
     *              + sigma2;
     */

    if( implicit )
      {
	const double denominator = ( 1 + uBar[0]*uBar[0]*k );
	const double denominator2 = denominator*denominator;

	ret[0] = -rho1 * ( u[1]*uBar[0]*uBar[0] / denominator + 2.0*u[0]*uBar[0]*uBar[1] / denominator2 )
	  + mu1 * u[0];
	ret[1] = rho2 * ( u[1]*uBar[0]*uBar[0] / denominator + 2.0*u[0]*uBar[0]*uBar[1] / denominator2 );

      }
    else
      {
	const double denominator = ( 1 + uBar[0]*uBar[0]*k );
	const double denominator2 = denominator*denominator;

	ret[0] = rho1 * ( -2.0 * uBar[0]*uBar[0]*uBar[1] / denominator2 );
	ret[1] = -rho2 * ( -2.0 *uBar[0]*uBar[0]*uBar[1] / denominator2 ) + sigma2;
      }
  }

private:
  const double a0, s0;
  const double var0;
  const double D1, D2;
  const double rho1, rho2, k, mu1, mu2, sigma1, sigma2;
};

template < class FunctionSpace >
class AllenCahnProblem : public TemporalProblemInterface < FunctionSpace >
{
 typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;
  using BaseType :: deltaT ;

  AllenCahnProblem( const TimeProviderType &timeProvider )
    : BaseType( timeProvider ),
      // initial values
      eps_ ( Dune::Fem::Parameter::getValue< double >( "allencahn.eps" ) )
  {}

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
    phi = RangeType(0);
    phi[ 0 ] = 0.1 * ( (double) rand() / RAND_MAX - 0.5 );
    phi[ 1 ] = 0.1 * ( (double) rand() / RAND_MAX - 0.5 );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& value) const
  {
    value = 0;
  }

  //! diffusion coefficient (default = Id)
  virtual void D( RangeType& D ) const
  {
    D = 0;

    // c1 diffusion coefficient
    D[ 0 ] = eps_;
    // c2 diffusion coefficient
    D[ 1 ] = eps_;
  }

  //! non linear term in reaction diffusion system
  virtual void nonLinearity( const RangeType &uBar,
			     const RangeType &u,
			     RangeType &ret,
			     const bool implicit ) const
  {
    /**
     *  full non-linearity is f(s) = ( s^3 - s ) / eps
     *
     *  this is linearized as
     *   f(s) approx f(s_0) + f'(s_0) ( s - s_0 )
     *  that is
     *   f(s) approx ( 3 s_0^2 - 1 ) * s - 2 * s0^3
     *
     *  here this is implemented as
     *  implicit:
     *    ret = u * ( -1 + 3 * uBar*uBar ) / eps
     *  explicit:
     *    ret = ( 2 uBar*uBar*uBar ) / eps
     */

    if( implicit )
      {
	for( unsigned int d = 0; d < ret.size(); ++d )
	  ret[ d ] = u[ d ] * ( - 1.0 + 3.0 * uBar[ d ] * uBar[ d ] );
	ret /= eps_;
      }
    else
      {
	for( unsigned int d = 0; d < ret.size(); ++d )
	  ret[ d ] = 2.0 * uBar[d]*uBar[d]*uBar[d];
	ret /= eps_;
      }
  }

private:
  const double eps_;
};

template < class FunctionSpace >
class BrusselatorProblem : public TemporalProblemInterface < FunctionSpace >
{
 typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;
  using BaseType :: deltaT ;

  BrusselatorProblem( const TimeProviderType &timeProvider )
    : BaseType( timeProvider ),
      alpha ( Dune::Fem::Parameter::getValue< double >( "Brusselator.alpha" ) ),
      beta ( Dune::Fem::Parameter::getValue< double >( "Brusselator.beta" ) ),
      gamma ( Dune::Fem::Parameter::getValue< double >( "Brusselator.gamma" ) ),
      r1 ( Dune::Fem::Parameter::getValue< double >( "Brusselator.r1" ) ),
      r2 ( Dune::Fem::Parameter::getValue< double >( "Brusselator.r2" ) ),
      D1 ( Dune::Fem::Parameter::getValue< double >( "Brusselator.D1" ) ),
      D2 ( Dune::Fem::Parameter::getValue< double >( "Brusselator.D2" ) )
  {}

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
    phi = RangeType(0);

#if TEST_CASE
    phi[ 0 ] = exp( -6.0 * D1 * time() ) * x[0] * x[1];
    phi[ 1 ] = exp( -6.0 * D2 * time() ) * x[0] * x[1];
#else
    // if( std::abs( x[0] ) < 0.2 )
    //   {
    // 	phi[ 0 ] = ( (double) rand() / RAND_MAX - 0.5 );
    // 	phi[ 1 ] = ( (double) rand() / RAND_MAX - 0.5 );
    //   }
    phi[ 0 ] = alpha + beta;
    phi[ 1 ] = beta / ((alpha+beta)*(alpha+beta));

    phi[ 0 ] += 0.1 * sin( 10 * x[0] ) * sin( 10 * x[1] ) * sin( 10 * x[2] );
    phi[ 1 ] += 0.1 * sin( 10 * x[0] ) * sin( 10 * x[1] ) * sin( 10 * x[2] );
#endif
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& value) const
  {
#if TEST_CASE
    const double u = exp( -6.0 * D1 * time() ) * x[0] * x[1];
    const double v = exp( -6.0 * D2 * time() ) * x[0] * x[1];

    value[ 0 ] = -( alpha * u * ( 1 - r1 * v*v ) + v * ( 1 - r2 * u ) );
    value[ 1 ] = -( v * ( beta + alpha * r1 * u * v ) + u * ( gamma + r2 * v ) );
#else
    value = 0;
#endif
  }

  //! diffusion coefficient (default = Id)
  virtual void D( RangeType& D ) const
  {
    D = 0;

    // c1 diffusion coefficient
    D[ 0 ] = D1;
    // c2 diffusion coefficient
    D[ 1 ] = D2;
  }

  //! non linear term in reaction diffusion system
  virtual void nonLinearity( const RangeType &uBar,
			     const RangeType &u,
			     RangeType &ret,
			     const bool implicit ) const
  {
    /**
     *  full non-linearity is f(u,v) = ( g * ( a - u + u^2 v ), g( b - u^2 v)
     *
     *  this is linearized as
     *   f(s) approx f(s_0) + f'(s_0) ( s - s_0 )
     *  that is
     *   f(s) approx ( 3 s_0^2 - 1 ) * s - 2 * s0^3
     *
     *  here this is implemented as
     *  implicit:
     *    ret = u * ( -1 + 3 * uBar*uBar ) / eps
     *  explicit:
     *    ret = ( 2 uBar*uBar*uBar ) / eps
     */

    if( implicit )
      {
	ret[0] = uBar[0]*uBar[0] * u[1] * gamma + u[0] * ( -1 + 2.0 * uBar[0]*uBar[1] ) * gamma;
	ret[1] = -uBar[0]*uBar[0] * u[1] * gamma + -2.0 * u[0] * uBar[0] * uBar[1] * gamma;
	ret *= -1.0; // now on rhs
      }
    else
      {
	// ret[ 0 ] = gamma * ( alpha - uBar[0] + uBar[0]*uBar[0] * uBar[1] );
	// ret[ 1 ] = gamma * ( beta - uBar[0]*uBar[0] * uBar[1] );
	ret[0] = ( -2.0 * uBar[0]*uBar[0]*uBar[1] + alpha ) * gamma;
	ret[1] = ( 2.0 *  uBar[0]*uBar[0]*uBar[1] + beta ) * gamma;
      }
  }

private:
  const double alpha, beta, gamma, r1, r2, D1, D2;
};

template < class FunctionSpace >
class GrayScottProblem : public TemporalProblemInterface < FunctionSpace >
{
 typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;
  using BaseType :: deltaT ;

  GrayScottProblem( const TimeProviderType &timeProvider )
    : BaseType( timeProvider ),
      // diffusion coefficients
      D1 ( Dune::Fem::Parameter::getValue< double >( "grayscott.D1" ) ),
      D2 ( Dune::Fem::Parameter::getValue< double >( "grayscott.D2" )),
      // initial values
      F ( Dune::Fem::Parameter::getValue< double >( "grayscott.F" ) ),
      k ( Dune::Fem::Parameter::getValue< double >( "grayscott.k" ))
  {}

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
    phi = RangeType(0);

    if( x[0] > 0.8 )
      {
	phi[ 0 ] = 0.5;
	phi[ 1 ] = 0.25;
      }
    else
      {
	phi[ 0 ] = 1.0;
	phi[ 1 ] = 0.0;
      }

    phi[ 0 ] += 0.1 * ( (double) rand() / RAND_MAX - 0.5 );
    phi[ 1 ] += 0.1 * ( (double) rand() / RAND_MAX - 0.5 );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& value) const
  {
    value = 0;
  }

  //! diffusion coefficient (default = Id)
  virtual void D( RangeType& D ) const
  {
    D = 0;

    // c1 diffusion coefficient
    D[ 0 ] = D1;
    // c2 diffusion coefficient
    D[ 1 ] = D2;
  }

  //! non linear term in reaction diffusion system
  virtual void nonLinearity( const RangeType &uBar,
			     const RangeType &u,
			     RangeType &ret,
			     const bool implicit ) const
  {
    if( implicit )
      {
	ret = 0;
      }
    else
      {
        ret[ 0 ] = - uBar[0]*uBar[1]*uBar[1] + F * ( 1 - uBar[0] );
        ret[ 1 ] = uBar[0]*uBar[1]*uBar[1] - ( F + k ) * uBar[1];
      }
  }

private:
  const double D1, D2, F, k;
};

#endif // #ifndef REACTIONDIFFUSION_HH
