function WriteRAMSIn(RamsData)


fid = fopen('rams.3.in', 'w' );   



fprintf( fid, 'rams.in\r\n\r\n');
fprintf( fid, '%6.2f\t 0\t %6.2f\t 1  \t ! freq rs  zs  direction\r\n', RamsData.freq, RamsData.zs  );

for ii = 1:length(RamsData.zr)
    fprintf( fid, '%6.2f\t', RamsData.zr(ii)  );
end;
fprintf( fid, '/receiver depths \r\n\r\n');

fprintf( fid, '%6.1f\t %2.4f\t %2.0f\t\t ! rmax dr  ndr\r\n', RamsData.rmax, RamsData.dr, RamsData.ndr );
fprintf( fid, '%6.1f\t %2.4f\t %2.0f\t %3.1f\t ! zmax	dz	ndz	zmplt\r\n', RamsData.zmax, RamsData.dz, RamsData.ndz, RamsData.zmplt );
fprintf( fid, '%6.1f\t %2.0f\t 0\t 0\t ! c0	np	irot	theta\r\n\r\n', RamsData.c0, RamsData.np);


for ii = 1:length(RamsData.bath)
    fprintf( fid, '%6.2f\t %6.3f    \r\n', RamsData.bath(ii,1), RamsData.bath(ii,2) );
end;
fprintf( fid, '-1 -1\r\n\r\n');


% if the file layers.txt exists, then profiles are constructed up to rmax
% if not, profiles are constructed up to the last hydrology profile
% the number of profiles corresponds to the cw size
% if bLayers has 1 column, then it is used for all profiles


nProf = size(RamsData.cw,2);
nzHydr = size(RamsData.cw,1);


for ii = 1:nProf
    
    if ii>1
        fprintf( fid, '%6.2f    \r\n', RamsData.drProf*(ii-1) );
    end;
    
    for jj = 1:nzHydr
        fprintf( fid, '%6.2f\t %6.3f    \r\n', RamsData.dzHydr*(jj-1), RamsData.cw(jj,ii) );
    end;
    fprintf( fid, '-1 -1\r\n\r\n');
    
    for jj = 1:5
        for kk = 1:size(RamsData.bLayers,1)
            if size(RamsData.bLayers,2) > 1
                fprintf( fid, '%6.2f\t %6.3f    \r\n', RamsData.bLayers(kk,ii), RamsData.bParams(kk,jj) );
            else
                fprintf( fid, '%6.2f\t %6.3f    \r\n', RamsData.bLayers(kk,1), RamsData.bParams(kk,jj) );
            end;
        end;
        
        fprintf( fid, '-1 -1\r\n\r\n');
    end;
    
    
end;


fclose( fid );