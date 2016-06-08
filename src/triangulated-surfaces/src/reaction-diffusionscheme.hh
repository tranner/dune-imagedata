#ifndef REACTIONDIFFUSION_FEMSCHEME_HH
#define REACTIONDIFFUSION_FEMSCHEME_HH

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
#include <dune/fem/solver/inverseoperators.hh>

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
#include "diffusionscheme.hh"

// Reaction-DiffusionScheme
//-------------------------
template < class ImplicitModel, class ExplicitModel >
struct ReactionDiffusionScheme 
  : public DiffusionScheme<ImplicitModel, ExplicitModel>
{
  typedef DiffusionScheme<ImplicitModel, ExplicitModel> BaseType;
  typedef typename BaseType::GridType GridType;
  typedef typename BaseType::GridPartType GridPartType;
  typedef typename BaseType::ModelType ImplicitModelType;
  typedef ExplicitModel ExplicitModelType;
  typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;

  ReactionDiffusionScheme( GridPartType &gridPart,
              ImplicitModelType& implicitModel,
              const ExplicitModelType& explicitModel )
    : BaseType(gridPart, implicitModel, explicitModel)
  {}

private:
  using BaseType::gridPart_;
  using BaseType::discreteSpace_;
  using BaseType::solution_;
  using BaseType::implicitModel_;
  using BaseType::rhs_;
  using BaseType::explicitModel_;
  using BaseType::linearOperator_;
  using BaseType::implicitOperator_;
  using BaseType::solverEps_;

  EllipticOperatorType explicitOperator_;
  std::vector< double > odeValues_;

  std::vector< double > solutionMassOld_;
  double surfaceAreaOld_;
  double interiorVolumeOld_;
};

#endif // end #if REACTIONDIFFUSION_FEMSCHEME_HH
