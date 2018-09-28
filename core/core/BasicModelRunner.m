classdef BasicModelRunner < IModelRunner
    
    properties
    end
    
    methods
        function  beforeLoadConfig (this, aModel)
            tic
            disp('beforeLoadConfig')
        end
        function  loadConfig (this, aModel)
            disp('loadConfig')
        end
        function  afterLoadConfig (this, aModel)
            disp('afterLoadConfig')
            elapsedTime = toc;
            fprintf('config loaded in %d ms\n', elapsedTime)
        end
        function  beforeRun (this, aModel)
            tic
            disp('beforeRun')
        end
        function  run (this, aModel)
            disp('run')
        end
        function  runUnix (this, aModel)
            disp('runUnix')
        end
        function  afterRun (this, aModel)
            disp('afterRun')
            elapsedTime = toc;
            fprintf('run executed in %d\n', elapsedTime)
        end
    end
end

