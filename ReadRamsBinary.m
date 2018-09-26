function [aField, dr, dz, varargout] = ReadRamsBinary(varargin)

pFolder = '';
if nargin > 0
    pFolder = varargin{1};
end;

dBounds = importdata([pFolder 'DomainBounds.Info'],':');
dr = dBounds.data(3);
dz = dBounds.data(6);
nr = dBounds.data(1);
nz = dBounds.data(4);

fid = fopen([pFolder 'TLrz']);
aField = fread(fid, [nz, nr],'float');
fclose(fid);

aField =  [ zeros(1,nr); aField];

if nargout>3
    fid = fopen([pFolder 'RePrz']);
    aFieldRe = fread(fid, [nz, nr],'float');
    fclose(fid);
    
    aFieldRe =  [ zeros(1,nr); aFieldRe];
    
    fid = fopen([pFolder 'ImPrz']);
    aFieldIm = fread(fid, [nz, nr],'float');
    fclose(fid);
    
    aFieldIm =  [ zeros(1,nr); aFieldIm];
    
    varargout{1} = aFieldRe + 1i*aFieldIm;
end;


