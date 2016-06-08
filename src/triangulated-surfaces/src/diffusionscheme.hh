#ifndef DIFFUSION_FEMSCHEME_HH
#define DIFFUSION_FEMSCHEME_HH

// iostream includes
#include <iostream>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// adaptation ...
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/common/adaptmanager.hh>

// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>

// include linear operators
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/solver/diagonalpreconditioner.hh>

#include <dune/fem/solver/istlsolver.hh>

// lagrange interpolation
#include <dune/fem/operator/lagrangeinterpolation.hh>

/*********************************************************/

// include norms
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// include parameter handling
#include <dune/fem/io/parameter.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

// local includes
#include "femscheme.hh"

// DiffusionScheme
//------------------

template < class ImplicitModel, class ExplicitModel >
struct DiffusionScheme : public FemScheme<ImplicitModel>
{
  typedef FemScheme<ImplicitModel> BaseType;
  typedef typename BaseType::GridType GridType;
  typedef typename BaseType::GridPartType GridPartType;
  typedef typename BaseType::ModelType ImplicitModelType;
  typedef ExplicitModel ExplicitModelType;
  typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
  DiffusionScheme( GridPartType &gridPart,
              const ImplicitModelType& implicitModel,
              const ExplicitModelType& explicitModel )
  : BaseType(gridPart, implicitModel),
    explicitModel_(explicitModel),
    explicitOperator_( explicitModel_, discreteSpace_ )
  {
  }

  void prepare()
  {
    // apply constraints, e.g. Dirichlet contraints, to the solution
    explicitOperator_.prepare( explicitModel_.dirichletBoundary(), solution_ );
    // apply explicit operator and also setup right hand side
    explicitOperator_( solution_, rhs_ );
    // apply constraints, e.g. Dirichlet contraints, to the result
    explicitOperator_.prepare( solution_, rhs_ );
  }

  void initialize ()
  {
     Dune::Fem::LagrangeInterpolation
          < typename ExplicitModelType::InitialFunctionType, DiscreteFunctionType > interpolation;
     interpolation( explicitModel_.initialFunction(), solution_ );
  }

  /** \brief compute interior volume of surface
   *
   *  computes the interior volume of surface using the formula
   *
   *  vol( Omega ) = 1/3 \int_\Gamma x \cdot \nu \, \mathrm{d} \sigma
   */
  double volume () const
  {
    typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    typedef typename IteratorType :: Entity EntityType;
    typedef typename EntityType :: Geometry GeometryType;
    typedef typename Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

    double vol = 0.0;

    const IteratorType end = discreteSpace_.end();
    for( IteratorType it = discreteSpace_.begin(); it != end; ++it )
      {
	const EntityType &entity = *it;
	const GeometryType &geometry = entity.geometry();

	// find normal
	Dune :: FieldVector < double, 3 > v0, v1, v2, normal;
	v0 = geometry.corner( 0 );
	v1 = geometry.corner( 1 );
	v2 = geometry.corner( 2 );

	v1 -= v0;
	v2 -= v0;

	normal[0] = v1[1]*v2[2]-v1[2]*v2[1];
	normal[1] = v1[2]*v2[0]-v1[0]*v2[2];
	normal[2] = v1[0]*v2[1]-v1[1]*v2[0];

	normal /= normal.two_norm();

	// perform quadrature
	const QuadratureType quadrature( entity, 5 );
	const size_t numQuadraturePoints = quadrature.nop();
	for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
	  {
	    // obtain quadrature point
	    const typename QuadratureType::CoordinateType &x = quadrature.point( pt );

	    Dune :: FieldVector < double, 3 > xGlobal = geometry.global( x );

	    // find integrand
	    const double integrand = xGlobal * normal;

	    vol += integrand * quadrature.weight( pt ) * geometry.integrationElement( x );
	  }
      }

    return vol / 3.0;
  }

protected:
  using BaseType::gridPart_;
  using BaseType::discreteSpace_;
  using BaseType::solution_;
  using BaseType::implicitModel_;
  using BaseType::rhs_;
  const ExplicitModelType &explicitModel_;
  typename BaseType::EllipticOperatorType explicitOperator_; // the operator for the rhs
};

#endif // end #if DIFFUSION_FEMSCHEME_HH
