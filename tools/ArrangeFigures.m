function ArrangeFigures(figs, rrin, ccin, offsetBottomin)
%ARRANGEFIGURES Выстраивает массив окон (figures) по экрану компьютера
%	figs	- массив указателей на объекты типа figures
%	rrcin	- количество строк в которые нужно построить окна
%	ccin	- количество столбцов в которые нужно построить окна
%	offsetBottomin - отступ снизу от низа экрана, чтобы учесть панель
%	задач. 
%
%Примеры вызовов :
%	ArrangeFigures(figs) - Выстроит массив окон figs в ряд по 
%	2 строки и 3 столбца. Отступ от низа экрана = 0;
%
%	ArrangeFigures(figs, rrin, ccin) - Выстроит массив окон figs в ряд по 
%	rrccin строк и ccin столбцов. Отступ от низа экрана = 0;
%
%	ArrangeFigures(figs, rrin, ccin, offsetBottomin) - Выстроит массив окон
%	figs в ряд по rrccin строк и ccin столбцов. Отступ от низа экрана = 
%	offsetBottomin пикселей
%
%Пример кода :
%	%Зададим массив объектов типа figures
%	hFig1=figure('Name','Figure1');
%	hFig2=figure('Name','Figure2');
%	hFig3=figure('Name','Figure3');
%	figs=[hFig1, hFig2, hFig3];
%	%Расположим эти объекты в 1 строку и три колонки на экране компьютера.
%	%Учтем, что панель задач перекрывает часть экрана компьютера, оценив ее 
%	%в 150 пикселей:
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