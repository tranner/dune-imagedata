#include <config.h>

// iostream includes
#include <iostream>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

// include deformable grid
#include <dune/grid/geometrygrid/grid.hh>
// include geometrty grid part
#include <dune/fem/gridpart/geogridpart.hh>
// include grid width
#include <dune/fem/misc/gridwidth.hh>
// include description of surface deformation
#include <dune/fem/function/blockvectorfunction.hh>
#include "deformation.hh"

// include header of adaptive scheme
#include "diffusion.hh"

#include "diffusionscheme.hh"
#include "reaction-diffusion.hh"
#include "reaction-diffusionmodel.hh"

// assemble-solve-estimate-mark-refine-IO-error-doitagain
template <class HGridType>
double algorithm ( HGridType &grid, int step )
{
  // create time provider
  typedef typename Dune::Fem::GridTimeProvider< HGridType > TimeProviderType;
  TimeProviderType timeProvider( grid );

  // choose problem
  typedef Dune::Fem::FunctionSpace< double, double, HGridType :: dimensionworld, 2 > FunctionSpaceType;
  typedef TemporalProblemInterface < FunctionSpaceType > ProblemInterfaceType;
  ProblemInterfaceType *problemPtr = 0;

  const std::string problemNames [] = { "KochMeinhardt", "AllenCahn", "Brusselator", "GrayScott" };
  const int problemNumber = Dune::Fem::Parameter :: getEnum("reactiondiffusion.problem", problemNames, 0);

  switch( problemNumber )
    {
    case 0:
      problemPtr = new KochMeinhardtProblem< FunctionSpaceType > ( timeProvider );
      break;

    case 1:
      problemPtr = new AllenCahnProblem< FunctionSpaceType > ( timeProvider );
      break;

    case 2:
      problemPtr = new BrusselatorProblem< FunctionSpaceType >( timeProvider );
      break;

    case 3:
      problemPtr = new GrayScottProblem< FunctionSpaceType >( timeProvider );
      break;

    default:
      std::cerr << "unrecognised problem name" << std::endl;
    }

  assert( problemPtr );
  ProblemInterfaceType &problem = *problemPtr;

  // choose type of discrete function
  typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimensionworld,
			       HGridType::dimensionworld > VectorFunctionSpaceType;
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType > HostGridPartType;
  HostGridPartType hostGridPart( grid );
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< VectorFunctionSpaceType, HostGridPartType, POLORDER > VertexSpaceType;
  VertexSpaceType vertexSpace( hostGridPart );
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< VertexSpaceType > VertexFunctionType;
  typedef DeformationDiscreteFunction< VertexFunctionType, TimeProviderType > DeformationType ;
  DeformationType deformation( vertexSpace, timeProvider );

#if !FIXED_MESH
  typedef Dune::Fem::GeoGridPart< DeformationType > GridPartType;
  GridPartType gridPart( deformation );
#else
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType, Dune::InteriorBorder_Partition > GridPartType;
  GridPartType gridPart( grid );
#endif

  const double hmax = Dune :: Fem :: GridWidth :: calcGridWidth( gridPart );
  if(  Dune :: Fem :: MPIManager ::rank() == 0 )
    std::cout << "mesh size h: " << hmax << std::endl;

  // type of the mathematical model used
  typedef ReactionDiffusionModel< FunctionSpaceType, GridPartType, DeformationType > ModelType;

  // implicit model for left hand side
  ModelType implicitModel( problem, gridPart, deformation, hmax, true );

  // explicit model for right hand side
  ModelType explicitModel( problem, gridPart, deformation, hmax, false );

  // create diffusion scheme
  typedef DiffusionScheme< ModelType, ModelType > SchemeType;
  SchemeType scheme( gridPart, implicitModel, explicitModel );

#if TEST_CASE
  typedef Dune::Fem::GridFunctionAdapter< ProblemInterfaceType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", problem, gridPart, 5 );
#endif

  //! input/output tuple and setup datawritter
#if TEST_CASE
  typedef Dune::tuple< const typename SchemeType::DiscreteFunctionType *, const GridExactSolutionType* > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(scheme.solution()), &gridExactSolution ) ; // tuple with pointers
#else
  typedef Dune::tuple< const typename SchemeType::DiscreteFunctionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(scheme.solution()) ) ; // tuple with pointers
#endif
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );

  const double endTime  = Dune::Fem::Parameter::getValue< double >( "reactiondiffusion.endtime", 100.0 );
  const double dtreducefactor = Dune::Fem::Parameter::getValue< double >("reactiondiffusion.reducetimestepfactor", 1 );
  double timeStep = Dune::Fem::Parameter::getValue< double >( "reactiondiffusion.timestep" );

  timeStep *= pow(dtreducefactor,step);

  // initialize with fixed time step
  timeProvider.init( timeStep ) ;

  // initialize scheme and output initial data
  scheme.initialize();
  // write initial solve
  dataOutput.write( timeProvider );
  timeProvider.next( timeStep );
  deformation.update();

  double nextMaxMinTime = 0;

  // time loop, increment with fixed time step
  for( ; timeProvider.time() < endTime; timeProvider.next( timeStep ) )
  {
    // find max min and print
    if( dataOutput.willWrite( timeProvider ) or timeProvider.time() > nextMaxMinTime )
    {
      double min0 = 1.0e6, min1 = 1.0e6;
      double max0 = -min0, max1 = -min1;
      int i = 0;
      for( auto dit = scheme.solution().dbegin();
	   dit != scheme.solution().dend(); ++dit )
	{
	  if( i++ % 2 == 0 )
	    {
	      min0 = std::min( *dit, min0 );
	      max0 = std::max( *dit, max0 );
	    }
	  else
	    {
	      min1 = std::min( *dit, min1 );
	      max1 = std::max( *dit, max1 );
	    }
	}
      std::cout << "time: " << timeProvider.time()
		<< " max0: " << max0
		<< " min0: " << min0
		<< " diff: " << max0 - min0 << std::endl;
      std::cout << "\t"
		<< " max1: " << max1
		<< " min1: " << min1
		<< " diff: " << max1 - min1 << std::endl;

      nextMaxMinTime += 1.0e-2;
    }

    // assemble explicit pare
    scheme.prepare();
    // move to new surface
    deformation.update();
    // solve once - but now we need to reassmble
    scheme.solve(true);

    dataOutput.write( timeProvider );
  }

  // output final solution
  dataOutput.write( timeProvider );

#if TEST_CASE
  // select norm for error computation
  typedef Dune::Fem::L2Norm< GridPartType > NormType;
  NormType norm( gridPart );
  return norm.distance( gridExactSolution, scheme.solution() );
#else
  return 0.0;
#endif
}

// main
// ----

int main ( int argc, char **argv )
try
{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize( argc, argv );

  // append overloaded parameters from the command line
  Dune::Fem::Parameter::append( argc, argv );

  // append possible given parameter files
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append( argv[ i ] );

  // append default parameter file
  Dune::Fem::Parameter::append( "../data/parameter" );

  // type of hierarchical grid
  //typedef Dune :: AlbertaGrid< 2 , 2 > GridType;
  typedef Dune :: GridSelector :: GridType  HGridType ;

  // create grid from DGF file
  const std::string gridkey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
  const std::string gridfile = Dune::Fem::Parameter::getValue< std::string >( gridkey );

  // the method rank and size from MPIManager are static
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro grid: " << gridfile << std::endl;

  // construct macro using the DGF Parser
  Dune::GridPtr< HGridType > gridPtr( gridfile );
  HGridType& grid = *gridPtr ;

  // do initial load balance
  grid.loadBalance();

  // setup EOC loop
  const int repeats = Dune::Fem::Parameter::getValue< int >( "reactiondiffusion.repeats", 0 );

  // initial grid refinement
  const int level = Dune::Fem::Parameter::getValue< int >( "reactiondiffusion.level" );

  // number of global refinements to bisect grid width
  const int refineStepsForHalf = Dune::DGFGridInfo< HGridType >::refineStepsForHalf();

  // refine grid
  Dune::Fem::GlobalRefine::apply( grid, level * refineStepsForHalf );

  // calculate first step
#if TEST_CASE
  double oldError = algorithm( grid, (repeats > 0) ? 0 : -1 );
  std::cout << "Error: " << oldError << std::endl;
  std::cout << std::endl;
#else
  algorithm( grid, (repeats > 0) ? 0 : -1 );
#endif

  for( int step = 1; step <= repeats; ++step )
    {
      // refine globally such that grid with is bisected
      // and all memory is adjusted correctly
      Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );

#if TEST_CASE
      const double newError = algorithm( grid, step );
      const double eoc = log( oldError / newError ) / M_LN2;
      if( Dune::Fem::MPIManager::rank() == 0 )
	{
	  std::cout << "Error: " << newError << std::endl;
	  std::cout << "EOC( " << step << " ) = " << eoc << std::endl;
	  std::cout << std::endl;
	}
      oldError = newError;
#else
      algorithm( grid, step );
#endif
    }

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
