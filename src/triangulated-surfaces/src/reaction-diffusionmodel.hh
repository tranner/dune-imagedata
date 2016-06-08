#ifndef REACTIONDIFFUSION_MODEL_HH
#define REACTIONDIFFUSION_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>

#include "temporalprobleminterface.hh"
#include "model.hh"

template< class FunctionSpace, class GridPart, class Deformation >
struct ReactionDiffusionModel : public EllipticModel<FunctionSpace,GridPart>
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

  static const unsigned int dimDomain = FunctionSpaceType :: dimDomain;
  static const unsigned int dimRange = FunctionSpaceType :: dimRange;

  typedef typename BaseType::ProblemType::DiffusionTensorType DiffusionTensorType;

  typedef typename DeformationType::VelocityType VelocityType;
  typedef typename DeformationType::NormalType NormalType;

  typedef TemporalProblemInterface< FunctionSpaceType > ProblemType ;

  typedef typename BaseType::ProblemType InitialFunctionType;

  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  //! constructor taking problem reference, time provider, 
  //! time step factor( either theta or -(1-theta) ), 
  //! flag for the right hand side 
  ReactionDiffusionModel( const ProblemType& problem, 
			  const GridPart &gridPart,
			  const DeformationType &deformation,
			  const double hmax,
			  const bool implicit )
    : BaseType(problem,gridPart),
      problem_( problem ),
      timeProvider_(problem.timeProvider()),
      deformation_( deformation ),
      hmax_( hmax ),
      implicit_( implicit ),
      timeStepFactor_( 0 ),
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

  //! this model does not have any mass term
  template< class Entity, class Point >
  void source ( const Entity &entity,
		const Point &x,
		const RangeType &value,
		RangeType &flux ) const
  {
    linSource( value, entity, x, value, flux );

    if( not implicit_ )
      {
	const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
	// evaluate right hand side
	RangeType rhs ;
	problem_.f( xGlobal, rhs );
	rhs  *= timeProvider_.deltaT();
	flux += rhs;
      }
  }

  //! this model does not have any mass term 
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

    // evaluate nonlinearity
    RangeType phi;
    problem_.nonLinearity( uBar, value, phi, implicit_ );
    // multiply by timestep and time step factor.
    phi *= timeProvider_.deltaT();
    flux += phi;
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
    RangeType DD;
    problem_.D( DD );
    for( int i =0; i < flux.rows; ++i )
      flux[ i ] *= DD[ i ];

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
	    for( int d = 0; d < dimRange; ++d )
	      {
		for( int i = 0; i < 3; ++i )
		  {
		    flux[ d ][ i ] += tangentialVelocity[ i ] * value[ d ];
		  }
	      }
	  }

	if( useStreamlineDiffusion_ )
	  {
	    // streamline diffusion
	    JacobianRangeType streamlineFlux(0);

	    // find velocity of surface
	    VelocityType velocity;
	    deformation_.velocity( x, entity, velocity );

	    for( int d = 0; d < dimRange; ++d )
	      {
		double dot2 = 0.0;
		for( int i = 0; i < 3; ++i )
		  dot2 += gradient[ d ][ i ] * velocity[ i ];

		for( int i = 0; i < 3; ++i )
		  streamlineFlux[ d ][ i ] = dot2 * velocity[ i ];
	      }

	    streamlineFlux *= hmax_*hmax_;
	    flux += streamlineFlux;
	  }
    }

    // multiply by timestep and time step factor.
    flux *= timeStepFactor_ * timeProvider_.deltaT();
  }

  //! return reference to Problem's time provider 
  const TimeProviderType & timeProvider() const 
  {
    return timeProvider_;
  }

  bool hasDirichletBoundary () const 
  {
    return false;
  }
  bool isDirichletPoint( const DomainType& x ) const 
  {
    return false;
  }
  const InitialFunctionType &initialFunction() const
  {
    return problem_;
  }

protected:
  const ProblemType &problem_;
  const TimeProviderType &timeProvider_; 
  const DeformationType &deformation_;
  const double hmax_;
  bool implicit_;
  double timeStepFactor_;
  const bool useStreamlineDiffusion_;
  const bool doSubtractTangentialVelocity_;
};

#endif // #ifndef REACTIONDIFFUSION_MODEL_HH
