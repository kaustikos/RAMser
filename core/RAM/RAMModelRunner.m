classdef RAMModelRunner < BasicModelRunner
    
    methods
        function  loadConfig (this, aModel)
            disp('RAM: loadConfig started')
            RamsData = LoadConfigRAMS(aModel.pFolder);
            WriteRAMSIn(RamsData);
        end
        function  run(this, aModel)
            disp('RAM: run started')
            system([aModel.rootFolder, '\bin\ram.exe'])
        end
        function  runUnix(this, aModel)
            disp('RAM: run started')
            system([aModel.rootFolder, '\bin\ram.bin'])
        end
        function  afterRun (this, aModel)
            afterRun@BasicModelRunner(this, aModel)            
            movefile('rams.3.in',[aModel.pFolder 'results/' 'rams.3.in']);
            movefile('DomainBounds.Info',[aModel.pFolder 'results/' 'DomainBounds.Info']);
            movefile('tl.nLine.Txt',[aModel.pFolder 'results/' 'tl.nLine.Txt']);
            movefile('TLrz',[aModel.pFolder 'results/' 'TLrz']);
            movefile('RePrz',[aModel.pFolder 'results/' 'RePrz']);
            movefile('ImPrz',[aModel.pFolder 'results/' 'ImPrz']);
        end
    end
end

