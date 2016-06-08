#ifndef DEFORMATION_HH
#define DEFORMATION_HH

#define CALIBRATION 0

#include <dune/common/exceptions.hh>
#include <dune/grid/geometrygrid/coordfunction.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

class FileNotFoundError : public Dune::IOError {};

// DeformationCoordFunction
// ------------------------

template < class DiscreteFunction, class TimeProvider >
class DeformationDiscreteFunction
  : public DiscreteFunction
{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef TimeProvider TimeProviderType;

public:
  typedef DiscreteFunction BaseType;

  typedef typename BaseType :: GridType GridType;
  typedef typename BaseType :: GridPartType GridPartType;
  typedef typename BaseType :: DiscreteFunctionSpaceType VertexFunctionSpaceType;
  typedef typename BaseType :: LocalFunctionType LocalFunctionType;

  typedef typename BaseType :: DomainType DomainType;
  typedef typename BaseType :: RangeType RangeType;

  typedef Dune::FieldVector< double, 3 > VelocityType;
  typedef Dune::FieldVector< double, 3 > NormalType;

  enum{ dimRange = 3 };

  DeformationDiscreteFunction ( const VertexFunctionSpaceType &vertexSpace,
				const TimeProviderType &timeProvider )
    : BaseType( "deformation", vertexSpace ),
      verticesNew_( "newDeformation", vertexSpace ),
      verticesOld_( "oldDeformation", vertexSpace ),
      timeProvider_( timeProvider ),
      nextFileTimeStep_( Dune::Fem::Parameter::getValue< double >( "imagedata.filetimestep", 0.1 ) ),
      nextFileTime_( 0.0 ),
      count_( Dune::Fem::Parameter::getValue< int >( "imagedata.initialfile", 0 ) )
#if TEST_CASE
    ,
      testEvolutionFunction_( timeProvider_ )
#endif
  {
#if TEST_CASE
    // initialise vertices
    Dune::Fem::LagrangeInterpolation
      < TestEvolutionFunction, DiscreteFunctionType > interpolation;
    interpolation( testEvolutionFunction_, (*this) );
    verticesOld_.assign( (*this) );
#else
#if CALIBRATION
    calibration();
#endif

#if !FIXED_MESH
    DGFToVertices( gridfile() );

    count_++;
    verticesOld_.assign( verticesNew_ );

    DGFToVertices( gridfile() );
    nextFileTime_ += nextFileTimeStep_;
    count_++;
#endif
#endif
    update();
  }

  DeformationDiscreteFunction ( const DeformationDiscreteFunction &other )
    : BaseType( other ),
      verticesOld_( other.verticesOld_ ),
      timeProvider_( other.timeProvider_ ),
      nextFileTime_( other.nextFileTime_ ),
      count_( other.count_ )
  {
    update();
  }

  void update ()
  {
#if !FIXED_MESH
  #if TEST_CASE
    verticesOld_.assign( (*this) );

    // interpolate test function at correct time
    Dune::Fem::LagrangeInterpolation
      < TestEvolutionFunction, DiscreteFunctionType > interpolation;
    interpolation( testEvolutionFunction_, (*this) );
    interpolation( testEvolutionFunction_, verticesNew_ );
  #else
    if( timeProvider_.time() > nextFileTime_ )
      {
	verticesOld_.assign( verticesNew_ );

	DGFToVertices( gridfile() );
	nextFileTime_ += nextFileTimeStep_;
	count_++;
      }

    // compute interpolation in time for mesh vertices
    {
      auto newIt = verticesNew_.dbegin();
      auto oldIt = verticesOld_.dbegin();
      auto dIt = this->dbegin();
      const auto dEnd = this->dend();

      const double l = ( timeProvider_.time() - ( nextFileTime_ - nextFileTimeStep_ ) ) / nextFileTimeStep_;

      for( ; dIt != dEnd; ++dIt, ++newIt, ++oldIt )
	{
	  *dIt = l * (*newIt) + (1-l) * (*oldIt);
	}
    }
  #endif
#endif
  }

#if TEST_CASE
  struct TestEvolutionFunction
    : public Dune :: AnalyticalCoordFunction< double, 3, 3, TestEvolutionFunction >
  {
    typedef Dune::Fem::FunctionSpace< double, double, 3, 3 > FunctionSpaceType;

    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    TestEvolutionFunction( const TimeProviderType &timeProvider )
      : timeProvider_( timeProvider )
    {}

    void evaluate( const DomainType &x, RangeType &y ) const
    {
      y = x;
#if !FIXED_MESH
      const double a = 1.0 + 0.25 * sin( time() );
      y[ 0 ] *= std::sqrt( a );

      const double theta = time() * M_PI;
      y[ 1 ] = x[ 1 ] * cos( theta ) - x[ 2 ] * sin( theta );
      y[ 2 ] = x[ 1 ] * sin( theta ) + x[ 2 ] * cos( theta );
#endif
    }

    double time() const
    {
      return timeProvider_.time();
    }

  private:
    const TimeProviderType &timeProvider_;
  };
#else

  std::string gridfile()
  {
    static const std::string dataBasename = Dune::Fem::Parameter::getValue< std::string >( "imagedata.databasename" );
    std::stringstream ss;
    ss << dataBasename << std::setw( 6 ) << std::setfill( '0' )
       << count_+1 << ".vertices";
    return ss.str();
  }

  void DGFToVertices( const std::string gridfile )
  {
    if( Dune::Fem::MPIManager::rank() == 0 )
      std::cout << "Loading vertex positions: " << gridfile << std::endl;

    auto dit = verticesNew_.dbegin();

    std::string line;
    std::ifstream file( gridfile );

    if( not file.is_open() )
      DUNE_THROW(FileNotFoundError, "File " << gridfile << " not found!");

    while( getline( file, line ) )
      {
	std :: stringstream ss;
	ss << line;

	for( int i = 0; i < 3; ++i )
	  {
	    double a;
	    ss >> a;
	    *dit = a;
	    ++dit;
	  }
      }
    file.close();
  }
#endif

  template< class HostEntity , class RangeVector >
  void evaluate ( const HostEntity &hostEntity, unsigned int corner,
                  RangeVector &y ) const
  {
    abort();
  }

  template <class RangeVector>
  void evaluate ( const typename GridType :: template Codim<0>::Entity &entity,
		  unsigned int corner,
                  RangeVector &y ) const
  {
    // find reference element
    typedef typename GridType::ctype  ctype;
    enum { dim = GridType::dimension };
    const Dune::ReferenceElement< ctype, dim > &refElement
      = Dune::ReferenceElements< ctype, dim >::general( entity.type() );

    // evaluate vertices
    LocalFunctionType localVertices = (*this).localFunction( entity );
    localVertices.evaluate( refElement.position( corner, dim ), y );
  }

  template< class Point, class Entity >
  void velocity ( const Point &pt, const Entity &entity,
		  VelocityType &velocity ) const
  {
#if !FIXED_MESH
    LocalFunctionType localVertices = verticesNew_.localFunction( entity.impl().hostEntity() );
    LocalFunctionType localVerticesOld = verticesOld_.localFunction( entity.impl().hostEntity() );

    // interpolate between two surfaces
    RangeType y, yold;
    localVertices.evaluate( pt, y );
    localVerticesOld.evaluate( pt, yold );

    velocity = ( y - yold );
    velocity /= timeProvider_.deltaT();
#else
    velocity = VelocityType(0);
#endif
  }

  template< class Point, class Entity >
  void normal ( const Point &pt, const Entity &entity,
		NormalType &normal ) const
  {
#if !FIXED_MESH
    LocalFunctionType localVertices = (*this).localFunction( entity.impl().hostEntity() );

    typedef typename GridType::ctype  ctype;
    enum { dim = GridType::dimension };

    const Dune::ReferenceElement< ctype, dim > &refElement
      = Dune::ReferenceElements< ctype, dim >::general( entity.type() );

    RangeType v0, v1, v2;
    localVertices.evaluate( refElement.position( 0, dim ), v0 );
    localVertices.evaluate( refElement.position( 1, dim ), v1 );
    v1 -= v0;
    localVertices.evaluate( refElement.position( 2, dim ), v2 );
    v2 -= v0;

    normal[0] = v1[1]*v2[2]-v1[2]*v2[1];
    normal[1] = v1[2]*v2[0]-v1[0]*v2[2];
    normal[2] = v1[0]*v2[1]-v1[1]*v2[0];

    normal /= normal.two_norm();
#else
    normal = NormalType(0);
#endif
  }

  const DiscreteFunctionType& vertices () const { return (*this); }
  DiscreteFunctionType& vertices () { return (*this); }

#if CALIBRATION
  void calibration()
  {
    const auto dend = this->dend();
    for( auto dit = this->dbegin(); dit != dend; ++dit )
      *dit = std::numeric_limits< double >::infinity();

    for( const auto& entity : this->space() )
      {
	const auto geometry = entity.geometry();
	const auto& lagrangePointSet = this->space().lagrangePointSet( entity );

	auto localVertices = this->localFunction( entity );

	// assume point based local dofs
        const int nop = lagrangePointSet.nop();
        int k = 0;
        for( int qp = 0; qp < nop; ++qp )
        {
          // if the first DoF for this point is already valid, continue
          if( localVertices[ k ] == std::numeric_limits< double >::infinity() )
          {
            // evaluate the function in the Lagrange point
            RangeType phi = geometry.global( coordinate( lagrangePointSet[ qp ] ) );

            // assign the appropriate values to the DoFs
            for( int i = 0; i < dimRange; ++i, ++k )
              localVertices[ k ] = phi[ i ];
          }
          else
            k += dimRange;
	}
      }

    std::ofstream file;
    file.open("calibration.txt");

    auto dit = (*this).dbegin();
    for( ; dit != (*this).dend(); ++dit )
      {
	file << *dit << "\n";
      }
    file << std::endl;

    assert(0);
  }
#endif

private:
  DiscreteFunctionType verticesNew_;
  DiscreteFunctionType verticesOld_;

  const TimeProviderType &timeProvider_;
  double nextFileTimeStep_;
  double nextFileTime_;
  int count_;

#if TEST_CASE
  TestEvolutionFunction testEvolutionFunction_;
#endif
};

#endif // #ifndef DEFORMATION_HH
