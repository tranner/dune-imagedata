#ifndef ELLIPTIC_HH
#define ELLIPTIC_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include "dirichletconstraints.hh"

// EllipticOperator
// ----------------

//! [Class for elliptic operator]
template< class DiscreteFunction, class Model >
struct EllipticOperator
: public virtual Dune::Fem::Operator< DiscreteFunction >
//! [Class for elliptic operator]
{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model            ModelType;

protected:
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType; 

  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  //! type of Dirichlet constraints 
  typedef Dune::DirichletConstraints< ModelType, DiscreteFunctionSpaceType > ConstraintsType;

public:
  //! contructor 
  EllipticOperator ( const ModelType &model, const DiscreteFunctionSpaceType &space )
  : model_( model )
  , constraints_( model, space )
  {}
      
  // prepare the solution vector 
  template <class Function>
  void prepare( const Function &func, DiscreteFunctionType &u ) 
  { 
    // set boundary values for solution 
    constraints()( func, u );
  }

  //! application operator 
  virtual void
  operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;

protected:
  const ModelType &model () const { return model_; }
  const ConstraintsType &constraints () const { return constraints_; }

private:
  ModelType model_;
  ConstraintsType constraints_;
};

// DifferentiableEllipticOperator
// ------------------------------

//! [Class for linearizable elliptic operator]
template< class JacobianOperator, class Model >
struct DifferentiableEllipticOperator
: public EllipticOperator< typename JacobianOperator::DomainFunctionType, Model >,
  public Dune::Fem::DifferentiableOperator< JacobianOperator >
//! [Class for linearizable elliptic operator]
{
  typedef EllipticOperator< typename JacobianOperator::DomainFunctionType, Model > BaseType;

  typedef JacobianOperator JacobianOperatorType;

  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename BaseType::ModelType ModelType;

protected:
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;

  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType; 

  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;

  typedef typename BaseType::QuadratureType QuadratureType;

public:
  //! contructor 
  DifferentiableEllipticOperator ( const ModelType &model, const DiscreteFunctionSpaceType &space, bool sw=true )
  : BaseType( model, space )
  {}
      
  //! method to setup the jacobian of the operator for storage in a matrix 
  void jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const;

protected:
  using BaseType::model;
  using BaseType::constraints;
};

// Implementation of EllipticOperator
// ----------------------------------

template< class DiscreteFunction, class Model >
void EllipticOperator< DiscreteFunction, Model >
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const 
{
  // clear destination 
  w.clear();

  // get discrete function space 
  const DiscreteFunctionSpaceType &dfSpace = w.space();

  // iterate over grid 
  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    // get entity (here element) 
    const EntityType &entity = *it;
    // get elements geometry 
    const GeometryType &geometry = entity.geometry();

    // get local representation of the discrete functions 
    const LocalFunctionType uLocal = u.localFunction( entity );
    LocalFunctionType wLocal = w.localFunction( entity );

    // obtain quadrature order 
    const int quadOrder = uLocal.order() + wLocal.order();

    { // element integral
      QuadratureType quadrature( entity, quadOrder );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        //! [Compute local contribution of operator]
        const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
        const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

        RangeType vu;
        uLocal.evaluate( quadrature[ pt ], vu );
        JacobianRangeType du;
        uLocal.jacobian( quadrature[ pt ], du );

        // compute mass contribution (studying linear case so linearizing around zero) 
        RangeType avu( 0 );
        model().source( entity, quadrature[ pt ], vu, avu );
        avu *= weight;
        // add to local functional wLocal.axpy( quadrature[ pt ], avu );

        JacobianRangeType adu( 0 );
        // apply diffusive flux 
        model().diffusiveFlux( entity, quadrature[ pt ], vu, du, adu );
        adu *= weight;

        // add to local function 
        wLocal.axpy( quadrature[ pt ], avu, adu );
        //! [Compute local contribution of operator]
      }
    }
  }

  // communicate data (in parallel runs)
  w.communicate();

  // apply constraints, e.g. Dirichlet contraints, to the result 
  constraints()( u, w );
}

// Implementation of DifferentiableEllipticOperator
// ------------------------------------------------

template< class JacobianOperator, class Model >
void DifferentiableEllipticOperator< JacobianOperator, Model >
  ::jacobian ( const DiscreteFunctionType &u, JacobianOperator &jOp ) const
{
  typedef typename JacobianOperator::LocalMatrixType LocalMatrixType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;

  const DiscreteFunctionSpaceType &dfSpace = u.space();

  Dune::Fem::DiagonalStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> stencil(dfSpace,dfSpace);
  jOp.reserve(stencil);
  jOp.clear();

  const int blockSize = dfSpace.localBlockSize; // is equal to 1 for scalar functions
  std::vector< typename LocalFunctionType::RangeType > phi( dfSpace.blockMapper().maxNumDofs()*blockSize );
  std::vector< typename LocalFunctionType::JacobianRangeType > dphi( dfSpace.blockMapper().maxNumDofs()*blockSize );

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    const LocalFunctionType uLocal = u.localFunction( entity );
    LocalMatrixType jLocal = jOp.localMatrix( entity, entity );

    const BasisFunctionSetType &baseSet = jLocal.domainBasisFunctionSet();
    const unsigned int numBasisFunctions = baseSet.size();
          
    QuadratureType quadrature( entity, 2*dfSpace.order() );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      //! [Assembling the local matrix]
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
      const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

      // evaluate all basis functions at given quadrature point 
      baseSet.evaluateAll( quadrature[ pt ], phi );

      // evaluate jacobians of all basis functions at given quadrature point
      baseSet.jacobianAll( quadrature[ pt ], dphi );

      // get value for linearization
      RangeType u0;
      JacobianRangeType jacU0;
      uLocal.evaluate( quadrature[ pt ], u0 );
      uLocal.jacobian( quadrature[ pt ], jacU0 );

      RangeType aphi( 0 );
      JacobianRangeType adphi( 0 );
      for( unsigned int localCol = 0; localCol < numBasisFunctions; ++localCol )
      {
        // if mass terms or right hand side is present 
        model().linSource( u0, entity, quadrature[ pt ], phi[ localCol ], aphi );

        // if gradient term is present 
        model().linDiffusiveFlux( u0, jacU0, entity, quadrature[ pt ], phi[ localCol ], dphi[ localCol ], adphi );

        // get column object and call axpy method 
        jLocal.column( localCol ).axpy( phi, dphi, aphi, adphi, weight );
      }
      //! [Assembling the local matrix]
    }
  } // end grid traversal 

  // apply constraints to matrix operator 
  constraints().applyToOperator( jOp );
  jOp.communicate();
}

#endif // #ifndef ELLIPTIC_HH
