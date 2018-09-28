classdef IModelRunner
    %IRUNMODEL Summary of this class goes here
    %   Detailed explanation goes here
       
    methods (Abstract)      
        result = beforeLoadConfig(arg, aModel)
        result = loadConfig(arg, aModel)
        result = afterLoadConfig(arg, aModel)
        result = beforeRun(arg, aModel)
        result = run(arg, aModel)
        result = runUnix(arg, aModel)
        result = afterRun(arg, aModel)
    end
end

