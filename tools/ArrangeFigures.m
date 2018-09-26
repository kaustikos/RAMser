function ArrangeFigures(figs, rrin, ccin, offsetBottomin)
%ARRANGEFIGURES ����������� ������ ���� (figures) �� ������ ����������
%	figs	- ������ ���������� �� ������� ���� figures
%	rrcin	- ���������� ����� � ������� ����� ��������� ����
%	ccin	- ���������� �������� � ������� ����� ��������� ����
%	offsetBottomin - ������ ����� �� ���� ������, ����� ������ ������
%	�����. 
%
%������� ������� :
%	ArrangeFigures(figs) - �������� ������ ���� figs � ��� �� 
%	2 ������ � 3 �������. ������ �� ���� ������ = 0;
%
%	ArrangeFigures(figs, rrin, ccin) - �������� ������ ���� figs � ��� �� 
%	rrccin ����� � ccin ��������. ������ �� ���� ������ = 0;
%
%	ArrangeFigures(figs, rrin, ccin, offsetBottomin) - �������� ������ ����
%	figs � ��� �� rrccin ����� � ccin ��������. ������ �� ���� ������ = 
%	offsetBottomin ��������
%
%������ ���� :
%	%������� ������ �������� ���� figures
%	hFig1=figure('Name','Figure1');
%	hFig2=figure('Name','Figure2');
%	hFig3=figure('Name','Figure3');
%	figs=[hFig1, hFig2, hFig3];
%	%���������� ��� ������� � 1 ������ � ��� ������� �� ������ ����������.
%	%�����, ��� ������ ����� ����������� ����� ������ ����������, ������ �� 
%	%� 150 ��������:
%	ArrangeFigures(figs, 1, 3, 150);

%% Initialization
	set(0,'Units','pixels') 
	scnsize = get(0,'ScreenSize');	
	counter = 0;
	rrinv = 0;
	if(nargin==1)
		RR=2;
		CC=3;
		offsetBottom=0;
	elseif(nargin==3)
		RR=rrin;
		CC=ccin;
		offsetBottom=0;
	elseif(nargin==4)
		RR=rrin;
		CC=ccin;
		offsetBottom=offsetBottomin;
	else
		error('unsupported numbers of input arguments');
	end
	
%% Arranging
	for rr=RR:-1:1
		rrinv=rrinv+1;
		for cc=1:CC
			counter=counter+1;
			if(counter<=length(figs))
				fig1=figs(counter);
				position = get(fig1,'Position');
				outerpos = get(fig1,'OuterPosition');
				borders = outerpos - position;
				edge = -borders(1)/2;			
				pos1 = [edge+(cc-1)*(1/CC)*scnsize(3),...%left
					(rr-1)*(1/RR)*(scnsize(4))+(1/RR)*offsetBottom*rrinv,...%bottom				
					(1/CC)*scnsize(3),...%width
					(1/RR)*(scnsize(4))-(1/RR)*offsetBottom-edge];%height
				set(fig1,'OuterPosition',pos1);
				figure(fig1);
			end
		end
	end	
end