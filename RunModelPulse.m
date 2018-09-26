function [t_r pt_r] = RunModelPulse(pFolder, zr, R, cmax)

set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial');

% LOAD PULSE DATA
xx = dlmread([pFolder 'ft.txt']);
dt = xx(2,1)-xx(1,1);
ft_s = xx(:,2);
% END PULSE DATA

nr = length(R);
ampl_tol = 0.05;

% LOAD RAMS
RamsData = LoadConfigRAMS(pFolder);

fRAMsMin = 2;
% END LOAD RAMS

nf = length(ft_s);
Tf = dt*(nf-1);
fs = 1/dt;

hF1 = figure();
plot( dt*(0:nf-1),ft_s,'linewidth',1.5 );
xlabel('t, s');
ylabel('p(t), s');
title('Source function (emitted signal)');
grid on;
ArrangeFigures(hF1,1,2,150);


NFFT = 4*(2^nextpow2(nf));

Spf = ifft(ft_s,NFFT);
f = fs/2*linspace(0,1,NFFT/2+1);
ampl_sp = max(abs(Spf));
ind_f_last = find( abs( Spf(1:NFFT/2+1) ) > ampl_sp*ampl_tol,1,'last');
disp('Upper bound of freq spectrum');
disp(f(ind_f_last) );
ind_f_first = find( abs( Spf(1:NFFT/2+1) ) > ampl_sp*ampl_tol,1,'first');
disp('Lower bound of freq spectrum');
disp(f(ind_f_first) );



hF1 = figure();
plot(f, abs(Spf(1:NFFT/2+1)),'linewidth',2,'color','black' );
grid on;
xlabel('f, Hz');
title('Frequency spectrum of the source function','linewidth',1.5);
xlim([0 f(ind_f_last)]);
ArrangeFigures(hF1,1,2,150);

Pomeg(1:nr,1:NFFT) = 0;

dorun = 1;
errs = [];

if dorun
    for ii = ind_f_first:ind_f_last
    
        
        if (abs(Spf(ii)) > ampl_sp*ampl_tol) && (f(ii)>fRAMsMin)
            disp(['Processing frequency f='  num2str(f(ii)) ]);
            
            
            
            %-----RUN RAMS
            
            RamsData.freq = f(ii);
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
            
            %         movefile('rams.3.in',[pFolder 'results/' 'rams.3.in']);
            %         movefile('DomainBounds.Info',[pFolder 'results/' 'DomainBounds.Info']);
            %         movefile('tl.nLine.Txt',[pFolder 'results/' 'tl.nLine.Txt']);
            %         movefile('TLrz',[pFolder 'results/' 'TLrz']);
            %         movefile('RePrz',[pFolder 'results/' 'RePrz']);
            %         movefile('ImPrz',[pFolder 'results/' 'ImPrz']);
            %
            %         [~, dr, dz, aFieldP] = ReadRamsBinary([pFolder 'results/']);
            
            
            [~, dr, dz, aFieldP] = ReadRamsBinary('');
            
            
            aFieldP = aFieldP/(4*pi);
            
            nr_out = size(aFieldP,2);
            nz_out = size(aFieldP,1);
            
            
            xx = dlmread('tl.nLine.Txt');
            if rem(ii,10) == 0
                close all;
                hh1 = figure();
                imagesc((0:nr_out-1)*dr/1000, (0:nz_out-1)*dz/1000 ,  trlo(aFieldP));
                caxis([-120  -60]);
                grid on;
                xlabel('x, km');
                ylabel('z, km');
                
                hh2 = figure;
                plot(xx(:,1)/1000, xx(:,3));
                xlim([0 R(end)/1000]);
                xlabel('x, km');
                ArrangeFigures([hh1 hh2],1,2,400);
                
            end;
            
            izr = zr/dz + 1;
            iR = R/dr;
            
            PHelm = interp2(aFieldP,iR,izr);
            
            %-----END RUN RAMS
            %
            %         %-----RUN MODES
            %         lambda = cmin/f(ii);
            %         dz = min(lambda/ppwlength,minmeshsize);
            %         z = 0:dz:H;
            %
            %         if nargin == 6
            %             [wnum, wmode] = ac_modesr(z,MP,f(ii) );
            %         else
            %             [wnum, wmode] = ac_modesr(z,MP,f(ii),0 );
            %         end;
            %
            %         nmod = length(wnum);
            %         PHelm(1:nr) = 0;
            %
            %         for mm = 1:nmod
            %             psiz = interp1(z,wmode(:,mm),zr);
            %             psizs = interp1(z,wmode(:,mm),zs);
            %
            %             PHelm(1:nr) = PHelm(1:nr) + (1i*exp(-1i*pi/4)/sqrt(8*pi) )*psiz*psizs*exp(1i*wnum(mm)*R)./sqrt(wnum(mm)*R);
            %
            %         end;
            %         %-----END RUN MODES
            
            
            
            
            
            Pomeg(1:nr,ii) = Spf(ii)*PHelm(1:nr);
            
            if sum(isnan(Pomeg(1:nr,ii)),1)>0
                Pomeg(1:nr,ii) = 0;
                errs = [errs; [ii f(ii) 1] ];
            end;
            
            npts = size(xx,1);
            if xx(npts,5) > 2*xx( floor(0.75*npts) ,5)
                Pomeg(1:nr,ii) = 0;
                errs = [errs; [ii f(ii) 2] ];
            end;
            
            
            %Pomeg(1:nr,NFFT + 2 - ii) = conj( Pomeg(1:nr,ii) );

        end;
    end;
    dlmwrite([pFolder 'sp_domain_rams.txt'],Pomeg,'delimiter','\t','precision',8);
    dlmwrite([pFolder 'errs.txt'],errs,'delimiter','\t','precision',8);


else
    Pomeg = dlmread([pFolder 'sp_domain_rams.txt'],'\t');
    
end;

figure;
hold all;
for jj=1:nr
    plot(f, abs(Pomeg(jj,1:NFFT/2+1) ),'linewidth',1.5 );
end;
grid on;
xlabel('f, Hz');
title('Frequency spectrum at the receivers');
xlim([0 f(ind_f_last)]);
%legend(['R = ' num2str(R(1)/1000,3) ' km'],['R = ' num2str(R(2)/1000,3) ' km'],['R = ' num2str(R(3)/1000,3) ' km']  );

%
% figure;
% hold all;
% for jj=1:nr
%     plot( abs(Pomeg(jj,:) ),'linewidth',1.5 );
% end;
% grid on;
% title('Frequency spectrum of the source function');

tau0s = R/cmax;

% multiply by exponentials

omeg = 2*pi*linspace(0,fs,NFFT+1);
for jj = 1:nr
    Pomeg(jj,1:NFFT) = Pomeg(jj,1:NFFT).*exp( -1i*tau0s(jj)*omeg(1:NFFT) );
end;

% conjugate in order to make ifft real

for ii = 2:ind_f_last
    k0 = 2*pi*f(ii)/RamsData.c0;
    %Pomeg(1:nr,ii) = exp(1i*k0*R.').*Pomeg(1:nr,ii);
    Pomeg(1:nr,NFFT + 2 - ii) = conj( Pomeg(1:nr,ii) );
end;

% perform ifft

pt_r = real(fft( Pomeg(:,1:NFFT).' ));
t_r(1:NFFT,1:nr) = repmat(dt*(0:NFFT-1).',1,nr);

for jj = 1:nr
    t_r(1:NFFT,jj) = t_r(1:NFFT,jj) + tau0s(jj);
end;



figure;
hold all;
for jj=1:nr
    plot(t_r(1:NFFT,jj), pt_r(1:NFFT,jj) ,'linewidth',1.5 );
end;
grid on;
xlabel('t, s');
title('Signals at the receivers');



%%
% 
% WriteRAMSIn(RamsData);
% 
% tic
% if isunix
%     if strcmp(RamsData.aModel,'RAM')
%         !./ram.bin
%     else
%         !./rams.bin
%     end;
%     
% else
%     if strcmp(RamsData.aModel,'RAM')
%         !ram.exe
%     else
%         !rams.exe
%     end;
% end;
% if ~exist([pFolder 'results/'],'dir')
%     mkdir(pFolder,'results')
% end;
% toc

delete('TLrz');
delete('RePrz');
delete('ImPrz');
%
% movefile('rams.3.in',[pFolder 'results/' 'rams.3.in']);
% movefile('DomainBounds.Info',[pFolder 'results/' 'DomainBounds.Info']);
% movefile('tl.nLine.Txt',[pFolder 'results/' 'tl.nLine.Txt']);
