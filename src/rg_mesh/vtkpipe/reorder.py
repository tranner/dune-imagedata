cfile = open('calibration.txt')
vertices = []
c = 0
for l in cfile.readlines():
    if c == 0:
        v = []

    try:
        z = l.split()[0]
        v.append( float(z) )
    except:
        1

    c = c + 1
    if c == 3:
        vertices.append( v )
        c = 0

tot = [0,0,0]
for v in vertices:
    for i in range( len(tot) ):
        tot[ i ] = tot[ i ] + v[ i ]
center = [ s / len(vertices) for s in tot ]

diam = 0
for i in range( len( vertices ) ):
    for j in range( i ):
        d2 = sum( (a-b)**2 for a,b in zip( vertices[i], vertices[j] ) )
        if d2 > diam*diam:
            diam = d2**0.5

def transform( x ):
    return [ ( x_ - c_ ) / diam for x_, c_ in zip( x, center ) ]

dfile = open('stack000001_rgmesh.dgf')
dlines = dfile.readlines()
oldvertices = []
reading = False
for l in dlines:
    if 'VERTEX' in l:
        reading = True
        continue
    if '#' in l:
        reading = False
        continue

    if reading:
        data = l.split()
        oldvertices.append( [ float(data[0]) , float(data[1]), float( data[2] ) ] )

order = [ oldvertices.index( l ) for l in vertices ]

rescale_file = open('stackbase_rgmesh.dgf','w')
vbool = False
for line in dlines:
    if 'VERTEX' in line:
        vbool = True
        rescale_file.write(line)
        continue
    if '#' in line:
        vbool = False

    if vbool:
        data = line.split()
        v = [ float(d) for d in data ]
        Tv = transform( v )
        rescale_file.write(' {0} {1} {2}\n'.format( Tv[0], Tv[1], Tv[2] ) )
    else:
        rescale_file.write(line)
rescale_file.close()


count = 1
for i in xrange(1,215):
    try:
        dfile = open('stack{0:06d}_rgmesh.dgf'.format(i) )
    except:
        print 'unable to open file', i
        continue

    try:
        oldvertices = []
        simplex = []
        reading = False
        sreading = False
        for l in dfile.readlines():
            if 'VERTEX' in l:
                reading = True
                continue
            if 'SIMPLEX' in l:
                sreading = True
                continue
            if '#' in l:
                reading = False
                sreading = False
                continue

            if reading:
                data = l.split()
                oldvertices.append( [ float(data[0]) , float(data[1]), float( data[2] ) ] )
            if sreading:
                data = l.split()
                simplex.append( [ int(data[0]), int(data[1]), int(data[2]) ] )

    except:
        print 'unable to read file', i

    print len(oldvertices)

    try:
        vfile = open( 'stack_rgmesh{0:06d}.vertices'.format( i ), 'w' )
        for o in order:
            v = transform( oldvertices[ o ] )
            vfile.write( '{0} {1} {2}\n'.format( v[0], v[1], v[2] ) )
        vfile.close()

        dfile = open( 'stack_rgmesh{0:06d}_scaled.dgf'.format( i ), 'w' )
        dfile.write( 'DGF\n\n')

        dfile.write( 'VERTEX\n' )
        for oldv in oldvertices:
            v = transform( oldv )
            dfile.write( '{0} {1} {2}\n'.format( v[0], v[1], v[2] ) )
        dfile.write('#\n\n')

        dfile.write( 'SIMPLEX\n' )
        for s in simplex:
            dfile.write( '{0} {1} {2}\n'.format( s[0], s[1], s[2] ) )
        dfile.write('#\n\n')

        dfile.close()
    except:
        print 'unable to write file', count
        continue
