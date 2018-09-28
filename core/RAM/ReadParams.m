function RAMs_params = ReadParams

%hydrology
hydro_files=filesearch('txt','Hydrology',1);

RAMs_params.hcp = zeros(0,1);
RAMs_params.hz  = zeros(0,1);
RAMs_params.hr  = zeros(0,1);

for i=1:length(hydro_files)
    
    Fname=char(hydro_files(i));
    fid=fopen(Fname,'r');
    x=fscanf(fid,'%d %d',[2 inf]);
    x=x';

    fclose(fid);
    
    RAMs_params.hz=[RAMs_params.hz x(:,1)];
    RAMs_params.hcp=[RAMs_params.hcp x(:,2)];
    
    if strfind(Fname,'\')~=0
        RAMs_params.hr=[RAMs_params.hr str2num(Fname(strfind(Fname,'\')+1:strfind(Fname,'.txt')-1))];
    else
        RAMs_params.hr=[RAMs_params.hr str2num(Fname(strfind(Fname,'/')+1:strfind(Fname,'.txt')-1))];
    end
    
end

 [RAMs_params.hr ix]=sort(RAMs_params.hr);
 RAMs_params.hz=RAMs_params.hz(:,ix);
 RAMs_params.hcp=RAMs_params.hcp(:,ix);

 %bathimetria
 RAMs_params.bath = zeros(0,1);
 
 if strfind(Fname,'\')~=0
        Fname='Bathymetria\bath.txt';
    else
        Fname='Bathymetria/bath.txt';
 end
   
    fid=fopen(Fname,'r');
    x=fscanf(fid,'%d %d',[2 inf]);
    RAMs_params.bath=x';
    
 %bottom layers
 
  layer_files=filesearch('txt','Layers',1);

RAMs_params.lr  = zeros(0,1); 
RAMs_params.lz  = zeros(0,1);
RAMs_params.lcp = zeros(0,1);
RAMs_params.lcs = zeros(0,1);
RAMs_params.lro = zeros(0,1);
RAMs_params.lap = zeros(0,1);
RAMs_params.las = zeros(0,1);


for i=1:length(layer_files)
    
    Fname=char(layer_files(i));
    fid=fopen(Fname,'r');
    x=fscanf(fid,'%d %d %d %d %d %d',[6 inf]);
    x=x';

    fclose(fid);
    RAMs_params.lz  = [RAMs_params.lz x(:,1)]; 
    RAMs_params.lcp = [RAMs_params.lcp x(:,2)]; 
    RAMs_params.lcs = [RAMs_params.lcs x(:,3)]; 
    RAMs_params.lro = [RAMs_params.lro x(:,4)]; 
    RAMs_params.lap = [RAMs_params.lap x(:,5)]; 
    RAMs_params.las = [RAMs_params.las x(:,6)]; 
    
    if strfind(Fname,'\')~=0
        RAMs_params.lr=[RAMs_params.lr str2num(Fname(strfind(Fname,'\')+1:strfind(Fname,'.txt')-1))];
    else
        RAMs_params.lr=[RAMs_params.lr str2num(Fname(strfind(Fname,'/')+1:strfind(Fname,'.txt')-1))];
    end
    
end

[RAMs_params.lr ix]=sort(RAMs_params.lr);
 
RAMs_params.lz  = RAMs_params.lz(:,ix);
RAMs_params.lcp = RAMs_params.lcp(:,ix);
RAMs_params.lcs = RAMs_params.lcs(:,ix);
RAMs_params.lro = RAMs_params.lro(:,ix);
RAMs_params.lap = RAMs_params.lap(:,ix);
RAMs_params.las = RAMs_params.las(:,ix);

%found the farest point 
size_temp=size(RAMs_params.bath);
R_max=RAMs_params.bath(size_temp(1,1),1);

%interpolation
step=500;
xi = 0:step:R_max; 

size_temp=size(RAMs_params.hcp);

%for Cp
x_temp  = zeros(0,1); 
for i=1:size_temp(:,1);   
yi = interp1(RAMs_params.hr,RAMs_params.hcp(i,:),xi); 
x_temp=[x_temp; yi];
end

RAMs_params.hcp=x_temp;

%for z
x_temp  = zeros(0,1); 
for i=1:size_temp(:,1);   
yi = interp1(RAMs_params.hr,RAMs_params.hz(i,:),xi); 
x_temp=[x_temp; yi];
end

RAMs_params.hz=x_temp;
RAMs_params.hr=xi;
