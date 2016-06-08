#ifndef DIFFUSION_MODEL_HH
#define DIFFUSION_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>

#include "temporalprobleminterface.hh"
#include "model.hh"

template< class FunctionSpace, class GridPart, class Deformation >
struct DiffusionModel : public EllipticModel<FunctionSpace,GridPart>
{
  typedef EllipticModel<FunctionSpace,GridPart> BaseType;
  typedef FunctionSpace FunctionSpaceType;
  typedef GridPart GridPartType;
  typedef Deformation DeformationType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename DeformationType::VelocityType VelocityType;
  typedef typename DeformationType::NormalType NormalType;

  typedef TemporalProblemInterface< FunctionSpaceType > ProblemType ;

  typedef typename BaseType::ProblemType InitialFunctionType;

  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  //! constructor taking problem reference, time provider,
  //! time step factor( either theta or -(1-theta) ),
  //! flag for the right hand side
  DiffusionModel( const ProblemType& problem,
		  const GridPart &gridPart,
		  const DeformationType &deformation,
		  const double hmax,
		  const bool implicit )
    : BaseType(problem,gridPart),
      timeProvider_(problem.timeProvider()),
      deformation_( deformation ),
      hmax_( hmax ),
      implicit_( implicit ),
      timeStepFactor_( 0 ),
      D_( Dune::Fem::Parameter::getValue< double >("diffusion.D", 1 ) ),
      useStreamlineDiffusion_( Dune::Fem::Parameter::getValue< bool >("diffusion.streamlinediffusion", 1 ) ),
      doSubtractTangentialVelocity_( Dune::Fem::Parameter::getValue< bool >("diffusion.subtracttangentialvelocity", 1 ) )
  {
    // get theta for theta scheme
    const double theta = Dune::Fem::Parameter::getValue< double >("diffusion.theta", 0.5 );
    if (implicit)
      timeStepFactor_ = theta ;
    else
      timeStepFactor_ = -( 1.0 - theta ) ;
  }

  template< class Entity, class Point >
  void source ( const Entity &entity,
                const Point &x,
                const RangeType &value,
                RangeType &flux ) const
  {
    linSource( value, entity, x, value, flux );
    // the explicit model should also evaluate the RHS
    if( !implicit_ )
    {
      const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
      // evaluate right hand side
      RangeType rhs ;
      problem_.f( xGlobal, rhs );
      rhs  *= timeProvider_.deltaT();
      flux += rhs ;
    }
  }

  template< class Entity, class Point >
  void linSource ( const RangeType& uBar,
                   const Entity &entity,
                   const Point &x,
                   const RangeType &value,
                   RangeType &flux ) const
  {
    const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
    RangeType m;
    problem_.m(xGlobal,m);
    for (unsigned int i=0;i<flux.size();++i)
      flux[i] = m[i]*value[i];
    flux *= timeStepFactor_ * timeProvider_.deltaT();
    // add term from time derivative
    flux += value;
  }

  //! return the diffusive flux
  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity,
                       const Point &x,
                       const RangeType &value,
                       const JacobianRangeType &gradient,
                       JacobianRangeType &flux ) const
  {
    linDiffusiveFlux( value, gradient, entity, x, value, gradient, flux );
  }
  //! return the diffusive flux
  template< class Entity, class Point >
  void linDiffusiveFlux ( const RangeType& uBar,
                          const JacobianRangeType &gradientBar,
                          const Entity &entity,
                          const Point &x,
                          const RangeType &value,
                          const JacobianRangeType &gradient,
                          JacobianRangeType &flux ) const
  {
    flux = gradient;
    flux *= D_;

    if( implicit_ )
      {
	if( doSubtractTangentialVelocity_ )
	  {
	    // find velocity of surface
	    VelocityType velocity;
	    deformation_.velocity( x, entity, velocity );

	    // find normal to surface
	    NormalType normal;
	    deformation_.normal( x, entity, normal );

	    // find tangential velocity of surface (to be removed)
	    double dot = 0.0;
	    for( int i = 0; i < 3; ++i )
	      {
		dot += velocity[ i ] * normal[ i ];
	      }

	    VelocityType tangentialVelocity;
	    for( int i = 0; i < 3; ++i )
	      {
		tangentialVelocity[ i ] = velocity[ i ] - dot * normal[ i ];
	      }

	    // remove tangential velocity
	    for( int i = 0; i < 3; ++i )
	      {
		flux[ 0 ][ i ] += tangentialVelocity[ i ] * value;
	      }
	  }

	if( useStreamlineDiffusion_ )
	  {
	    // streamline diffusion
	    JacobianRangeType streamlineFlux(0);

	    // find velocity of surface
	    VelocityType velocity;
	    deformation_.velocity( x, entity, velocity );

	    double dot2 = 0.0;
	    for( int i = 0; i < 3; ++i )
	      dot2 += gradient[ 0 ][ i ] * velocity[ i ];

	    for( int i = 0; i < 3; ++i )
	      streamlineFlux[ 0 ][ i ] = dot2 * velocity[ i ];

	    streamlineFlux *= hmax_*hmax_;
	    flux += streamlineFlux;
	  }
    }

    // multiply by timestep and time step factor.
    flux *= timeStepFactor_ * timeProvider_.deltaT();
  }

  //! exact some methods from the problem class
  bool hasDirichletBoundary () const
  {
    return BaseType::hasDirichletBoundary() ;
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  bool isDirichletPoint( const DomainType& x ) const
  {
    return BaseType::isDirichletPoint(x) ;
  }

  template< class Entity, class Point >
  void g( const RangeType& uBar,
          const Entity &entity,
          const Point &x,
          RangeType &u ) const
  {
    BaseType::g(uBar,entity,x,u);
  }

  // return Fem :: Function for Dirichlet boundary values
  typename BaseType::DirichletBoundaryType dirichletBoundary( ) const
  {
    return BaseType::dirichletBoundary();
  }

  //! return reference to Problem's time provider
  const TimeProviderType & timeProvider() const
  {
    return timeProvider_;
  }

  const InitialFunctionType &initialFunction() const
  {
    return problem_;
  }

protected:
  using BaseType::problem_;
  const TimeProviderType &timeProvider_;
  const DeformationType &deformation_;
  const double hmax_;
  bool implicit_;
  double timeStepFactor_;
  double D_;
  const bool useStreamlineDiffusion_;
  const bool doSubtractTangentialVelocity_;
};
#endif // #ifndef DIFFUSION_MODEL_HH
