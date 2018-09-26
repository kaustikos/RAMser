function outputDir = GetDirs(inputDirs, prefix, sortOrder)
%GETDIRS Получить список папок с возможностью сортировки по значению

%% Initialization
	if(nargin==1)
		sortOrder='ASC';
		prefix='';
		sortByValue=0;
	elseif(nargin==2)
		sortOrder='ASC';
		sortByValue=1;
	elseif(nargin==3)
		sortByValue=1;
	else
		error('unsupported numbers of input arguments');
	end
	
%% Execute
	outputDir = GetDirsOrFiles(inputDirs, 1, sortByValue, prefix, sortOrder);		
end

