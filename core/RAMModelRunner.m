classdef RAMModelRunner < BasicModelRunner
    
    methods
        function  beforeLoadConfig (this, aModel)
            beforeLoadConfig@BasicModelRunner(this, aModel)
            disp('RAMS!!!!!!')
        end       
    end
end

