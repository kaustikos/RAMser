function outputDir = GetDirsOrFiles(inputDirs, isdir, sortByValue, prefix, sortOrder)
%GETDIRSORFILES Получить список файлов или папок с возможностью сортровки по значению
% Пример использования GetDirs('c:\test\','DESC','F');
% будут считаны все подпапки из папки 'c:\test\'. У этих папок будут удалены
% все вхождения символа 'F' и затем будет отсортированный список по
% убыванию.
% Если необходимо просто считать папки без системных папок '.' '..'
% используйте функцию без 2го и 3го параметров.
% Первым параметром может быть как строка с указанием папки так и массив
% структур типа dir

%% Iniitialization
	if(strcmp(sortOrder,'ASC')||strcmp(sortOrder,'ascend'))
		sortOrder='ascend';
	elseif(strcmp(sortOrder,'DESC')||strcmp(sortOrder,'descend'))
		sortOrder='descend';
	elseif(isempty(sortOrder))
		sortOrder='ascend';
	else
		error(['second argument can has only ''ASC'', ''ascend'',''DESC''',...
			' or ''descend'' values. Or let it empty']);
	end

	if(ischar(inputDirs))
		inputDirs=dir(inputDirs);
	elseif(~isstruct(inputDirs))
		error('first input argument must be a string or a struct');
	end

	if(~isempty(prefix) && ~ischar(prefix))
		error('third argument must be a string');
	end

%% Get values
	values={};
	indexes=[];
	if(isdir)
		for kk = 1:length(inputDirs);
			if(inputDirs(kk).isdir && ~strcmp(inputDirs(kk).name,'.')...
					&& ~strcmp(inputDirs(kk).name,'..'))
				values{end+1}=inputDirs(kk).name;
				indexes=[indexes, kk];
			end
		end
	else
		for kk = 1:length(inputDirs);
			if(~inputDirs(kk).isdir && ~strcmp(inputDirs(kk).name,'.')...
					&& ~strcmp(inputDirs(kk).name,'..'))
				values{end+1}=inputDirs(kk).name;
				indexes=[indexes, kk];
			end
		end
	end

%% Sorting
	outputDir=[];
	if(~sortByValue)
		for kk = 1:length(indexes)
			outputDir = [outputDir inputDirs(indexes(kk))];
		end
	else
		[val, ind] =SortByValue(values,prefix,sortOrder);
		for kk = 1:length(values)
			outputDir = [outputDir inputDirs(indexes(ind(kk)))];
		end
	end
end

