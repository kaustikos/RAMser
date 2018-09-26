function names = ParseParams(filename)
filetext = fileread(filename);
%exp = '(?<key>\w+)\s*=\s*(?<value>[\d\w\'']+)[ ]*;[.]*[\n]{1}';
%doesn't account comment with ;... f.e. %r=0; - is also
%parsable parameter
exp = '(?<key>\w+)\s*=\s*(?<value>-*[\d\w\''[]., ]+)[ ]*;';
names = regexp(filetext, exp, 'names');
% for ii = 1:size(names, 2)
%     if (isprop(this, names(ii).key))
%         if (isnumeric(getfield(this, names(ii).key)))
%             if (~isempty(strfind(names(ii).value, '[')))
%                 setfield(this, names(ii).key, this.GetMassiveFromText(names(ii).value));
%             else
%                 setfield(this, names(ii).key, str2double(names(ii).value));
%             end
%         else
%             setfield(this, names(ii).key, strrep(names(ii).value, '''', ''));
%         end
%     end
% end
% end