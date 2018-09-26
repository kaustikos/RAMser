function RunModel(pFolder)

RamsData = LoadConfigRAMS(pFolder);
WriteRAMSIn(RamsData);

tic
if isunix
  if strcmp(RamsData.aModel,'RAM')
      !./ram.bin
  else
      !./rams.bin
  end;

else
  if strcmp(RamsData.aModel,'RAM')
      !ram.exe
  else
      !rams.exe
  end;
end;
if ~exist([pFolder 'results/'],'dir')
    mkdir(pFolder,'results')
end;
toc

movefile('rams.3.in',[pFolder 'results/' 'rams.3.in']);
movefile('DomainBounds.Info',[pFolder 'results/' 'DomainBounds.Info']);
movefile('tl.nLine.Txt',[pFolder 'results/' 'tl.nLine.Txt']);
movefile('TLrz',[pFolder 'results/' 'TLrz']);
movefile('RePrz',[pFolder 'results/' 'RePrz']);
movefile('ImPrz',[pFolder 'results/' 'ImPrz']);
