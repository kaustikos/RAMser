classdef Launcher
    %LAUNCHER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        aModelRunner
        aModel
    end
    
    methods
        function obj = Launcher(iModelRunner, iModel)
            if(~isa(iModelRunner, 'IModelRunner'))
                error('first arg must be of type IModelRunner')
            end
            if(~isa(iModel, 'IModel'))
                error('second arg must be of type IModel')
            end
            obj.aModelRunner = iModelRunner;
            obj.aModel = iModel;
        end
        
        function outputArg = run(this)
            this.aModelRunner.beforeLoadConfig(this.aModel);
            this.aModelRunner.loadConfig(this.aModel);
            this.aModelRunner.afterLoadConfig(this.aModel);
            this.aModelRunner.beforeRun(this.aModel);
            if isunix
                this.aModelRunner.runUnix(this.aModel);
            else
                this.aModelRunner.run(this.aModel);
            end
            this.aModelRunner.afterRun(this.aModel);            
        end
    end
end

