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
// include description of surface deformation
#include <dune/fem/function/blockvectorfunction.hh>
#include "deformation.hh"

// include header of adaptive scheme
#include "diffusion.hh"

#include "diffusionmodel.hh"
#include "diffusionscheme.hh"

template< class GridPart >
double computeMeshSize( const GridPart &gridPart )
{
  typedef typename GridPart::template Codim< 0 >::IteratorType IteratorType;
  typedef typename GridPart::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPart::template Codim< 0 >::GeometryType GeometryType;

  double hmax = 0;

  // iterate over the interior of the grid
  const IteratorType end = gridPart.template end< 0 >();
  for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
    {
      // obtain a reference to the entity
      const EntityType &entity = *it;
      const GeometryType geometry = entity.geometry();

      for( int i = 0; i < geometry.corners(); ++i )
	{
	  for( int j = 0; j < i; ++j )
	    {
	      hmax = std::max( hmax, ( geometry.corner(i) - geometry.corner(j) ).two_norm() );
	    }
	}
    }

  return hmax;
}


// assemble-solve-estimate-mark-refine-IO-error-doitagain
template <class HGridType>
double algorithm ( HGridType &grid, int step )
{
  typedef Dune::Fem::FunctionSpace< double, double, HGridType :: dimensionworld, 1 > FunctionSpaceType;

  // create time provider
  typedef typename Dune::Fem::GridTimeProvider< HGridType > TimeProviderType;
  TimeProviderType timeProvider( grid );

  // choose type of discrete function
  typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimensionworld,
			       HGridType::dimensionworld > VectorFunctionSpaceType;
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType, Dune::InteriorBorder_Partition > HostGridPartType;

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

  // initialise hmax
  double hmax = computeMeshSize( gridPart );
  
  // type of the mathematical model used
  typedef DiffusionModel< FunctionSpaceType, GridPartType, DeformationType > ModelType;
  typedef typename ModelType :: ProblemType ProblemType;

  // choose problem
  ProblemType* problemPtr = 0;
  const std::string problemNames [] = { "diffusion", "frap", "ligand" };
  const int problemNumber = Dune :: Fem :: Parameter :: getEnum( "diffusion.problem", problemNames );
  switch( problemNumber )
    {
    case 0:
      problemPtr = new DiffusionProblem< FunctionSpaceType > ( timeProvider );
      break;
    case 1:
      problemPtr = new FrapProblem< FunctionSpaceType > ( timeProvider );
      break;
    case 2:
      problemPtr = new LigandProblem< FunctionSpaceType > ( timeProvider );
      break;

    default:
      std::cerr << "unrecognised problem name" << std::endl;
    }

  // recover problem
  assert( problemPtr );
  ProblemType& problem = *problemPtr;

  // implicit model for left hand side
  ModelType implicitModel( problem, gridPart, deformation, hmax, true );

  // explicit model for right hand side
  ModelType explicitModel( problem, gridPart, deformation, hmax, false );

  // create diffusion scheme
  typedef DiffusionScheme< ModelType, ModelType > SchemeType;
  SchemeType scheme( gridPart, implicitModel, explicitModel );

  typedef Dune::Fem::GridFunctionAdapter< ProblemType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", problem, gridPart, 5 );
  //! input/output tuple and setup datawritter
  typedef Dune::tuple< typename SchemeType::DiscreteFunctionType *, GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(scheme.solution()), &gridExactSolution) ; // tuple with pointers
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );

  const double endTime  = Dune::Fem::Parameter::getValue< double >( "diffusion.endtime", 2.0 );
  const double dtreducefactor = Dune::Fem::Parameter::getValue< double >("diffusion.reducetimestepfactor", 1 );
  double timeStep = Dune::Fem::Parameter::getValue< double >( "diffusion.timestep", 0.125 );

  // open fire for concentration data
  std::ofstream concentrationFile;
  {
    std::string name = Dune :: Fem ::Parameter :: commonOutputPath() + "/";
    // add prefix for data file
    // name += Dune :: Fem :: Parameter :: getValue< std::string >( "fem.prefix" );
    name += "concentration.txt";
    std::cout << "writing concentration data to " << name << std::endl;
    concentrationFile.open( name );
  }

  if( step > 0 )
    timeStep *= pow(dtreducefactor,step);

  // initialize with fixed time step
  timeProvider.init( timeStep ) ;

  // initialize scheme and output initial data
  scheme.initialize();
  // write initial solve
  dataOutput.write( timeProvider );
  concentrationFile << timeProvider.time()
		    << " " << scheme.regionConcentration()
		    << " " << scheme.regionArea() << std::endl;

  // increment time step
  timeProvider.next( timeStep );
  deformation.update();

  // time loop, increment with fixed time step
  for( ; timeProvider.time() < endTime; timeProvider.next( timeStep ) )
  {
    // assemble explicit pare
    scheme.prepare();
    // move to new surface
    deformation.update();
    // solve once - but now we need to reassmble
    scheme.solve(true);

    dataOutput.write( timeProvider );
    concentrationFile << timeProvider.time()
		      << " " << scheme.regionConcentration()
		      << " " << scheme.regionArea() << std::endl;
  }

  // output final solution
  dataOutput.write( timeProvider );

  // close concentration file
  concentrationFile.close();

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
  const int repeats = Dune::Fem::Parameter::getValue< int >( "diffusion.repeats", 0 );

  // initial grid refinement
  const int level = Dune::Fem::Parameter::getValue< int >( "diffusion.level" );

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

