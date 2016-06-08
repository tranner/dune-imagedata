% DGFTRISURF % DGFQUIVER creates a dgf-file containing a plot similar to
%    Matlab's trisurf which can be called by a dgf-visualization tool.
%
%    DGFTRISURF(TRI,X,Y,Z,VALS,VARNAME,FILENAME)
%    See Matlab documentation acc. to TRISURF(TRI,X,Y,Z,C).
%
%    DGFTRISURF(TRI,X,Y,Z,VARNAME,FILENAME)
%    See Matlab documentation acc. to TRISURF(TRI,X,Y,Z).
%
%    The (optional) argument VALS should contain the scalar data, which can
%    be given on the vertices or on the cells/triangles of the triangu-
%    lation.  The type is determined automatically.If VALS is not 
%    specified, Z is used as point data.  The argument VARNAME should 
%    contain the description of the visualized variable which is required
%    by the dgf file format. FILENAME is the name of the .dgf-file to store
%    the data.  If the string does not end with '.dgf', this extension is 
%    attached.
%
% Example
%    [X,Y] = meshgrid(-2:0.25:2,-1:0.2:1);
%    Z = X.* exp(-X.^2 - Y.^2);
%    TRI = delaunay(X,Y);
%    dgftrisurf(TRI,X,Y,Z,'pressure','dgftrisurf.dgf')
%
% See also trisurf, vtkquiver

function dgftrisurf(tri, x, y, z, argin5, argin6, argin7)

% INITIALIZATION
if nargin == 6
  vals     = z;
  varname  = argin5;
  filename = argin6;
elseif nargin == 7
  vals     = argin5;
  varname  = argin6;
  filename = argin7;
else
  error('Wrong number of input arguments.')
end

% ASSERTIONS
assert(ischar(varname) && ischar(filename))
numC = size(tri, 1); % number of cells
numP = length(x(:)); % number of points
assert(numP == length(y(:)) && numP == length(z(:)))

% DETERMINE DATA TYPE
if length(vals(:)) == numP
  datatype = 'pointdata';
elseif length(vals(:)) == numC
  datatype = 'celldata';
else
  error('Input argument VALS has wrong dimensions.')
end

% OPEN FILE
if ~strcmp(filename(end-3:end), '.dgf') % append file extension if not specified yet
  filename = [filename '.dgf'];
end

file = fopen(filename, 'wt');

% HEADER
fprintf(file, 'DGF\n\n');

% VERTEX
fprintf(file, 'VERTEX\n');
for kV = 1 : numP
  fprintf(file, ' %.3e %.3e %.3e\n', x(kV), y(kV), z(kV));
end
fprintf(file, '#\n\n');

% SIMPLEX
fprintf(file, 'SIMPLEX\n');
for kC = 1 : numC
  fprintf(file,'           %d %d %d\n', tri(kC, 1) - 1, tri(kC, 2) - 1, tri(kC, 3) - 1);
end
fprintf(file, '#\n\n');

fclose(file);

return

end
