% ML_TimeSeriesUS - Plot US time series data with NBER Recessions
% 
% ML_TimeSeriesUS(ZZ,Dates,Titolo,legenda,asseX,sizeXY,ng,colore,ny,mq,ypace)
% 
%      ZZ = Variables to be plotted
%   Dates = Dates upon which ZZ has to be plotted (lenght(ZZ)==length(Dates)
%  Titolo = Title of the graph (not optional but can be empty)
% legenda = Entries for legend (not optional but can be empty)
%   asseX = xlabel (not optional but can be empty)
%  sizeXY = size of axis (optional)
%      ng = if ng=1 ==> no grid (optional)
%  colore = colors for the graph (optional)
%      ny = interval (in years) for grid on x-axis (optional)
%      mq = set in which month the grid of the x-axis is put (optional)
%   ypace = Pace for y-ticks (optional)

% Matteo Luciani (matteoluciani@yahoo.it)

function pl=ML_TimeSeriesUSbw(ZZ,Dates,legenda,sizeXY,LS,ny,mq,ypace)

%%% ===================================== %%%
%%%  Check inputs and set default values  %%%
%%% ===================================== %%%
if ~isempty(legenda); lg=1; else; lg=0; end
try isnan(sizeXY); catch; sizeXY=[]; end

try isnumeric(ny); if isempty(ny); ny=5; end;  catch; ny=5; end
try isnan(mq); if isempty(mq); mq=3; end; catch; mq=3; end
try isempty(ypace); ndec=['%.' num2str(ML_ndecimal(ypace)) 'f'];
    if isempty(ypace); ypace=1; ndec='%.0f'; end
catch, ypace=1; ndec='%.0f';
end

if isempty(LS); LS={'k-','k--','k:','k-.'}; end
[T,N]=size(ZZ);    

cww=[.9,.9,.9];                                                      % color for NBER recession


if isempty(sizeXY)                                                          % if size of figure not provided
    xlim=Dates([1 end])'; ylim=[ML_min(ZZ) ML_min(ZZ,2)];                   % ------------------------------
else; xlim=sizeXY(1:2); ylim=sizeXY(3:4);                                    % size of figure provided by User
end

if sum(sign(ylim))==0   % ------------------------------------------------- % Set-up yticks
    t1=0:-ypace:ML_round(ylim(1),0); t2=0:ypace:ML_round(ylim(2),0);
    jjj=0;
    while ML_round(ylim(1),jjj)==0&&ML_round(ylim(2),jjj)==0; jjj=jjj+1;
        t1=0:-ypace:ML_round(ylim(1),jjj); t2=0:ypace:ML_round(ylim(2),jjj);
    end
    ytick=sort([t1(2:end) t2]); 
elseif sum(sign(ylim))>0
    ytick=ML_round(ylim(1),0,1):ypace:ML_round(ylim(2),0);       
else; ytick=ML_round(ylim(1),1):ypace:ML_round(ylim(2),0);
end                     % ------------------------------------------------- %

temp=datevec(Dates); TICK=find(temp(:,2)==mq);                              % Tick at each mq of the year
TICK2=TICK(mod(temp(TICK,1),ny)==0);                                        % grid every ny years
xtl=cellstr(repmat(' ',length(TICK),1));                                    % define the stile of the ticklabel
xtl(mod(temp(TICK,1),ny)==0)=cellstr(datestr(Dates(TICK2),'yyyy'));         % ---------------------------------

load NBER_Recessions; NBER2=NBER(:,2)*ylim;

        %%% ======================= %%%
        %%%         Graphing        %%%
        %%% ======================= %%%
        
axes('Parent',figure,'FontSize',12); ML_FigureSize,hold on;
ha=area(NBER(:,1),[NBER2(:,1) NBER2(:,2)-NBER2(:,1)],'linestyle','none');  	% NBER Recessions
set(ha(1), 'FaceColor', 'none'); set(ha(2), 'FaceColor', cww)               % ---------------
axis([xlim ylim])                                                           % define size of the figure
set(gca,'Xtick', Dates(TICK),'Xticklabel',xtl);                             % set xtick
set(gca,'Ytick',ytick,'Yticklabel',num2str(ytick',ndec));                   % set ytick
plot(Dates,zeros(T,1),'k','linewidth',.5);                                  % line at zero
for nn=1:N  % ------------------------------------------------------------- % plot each variable separately
    pl(nn)=plot(Dates,ZZ(:,nn),LS{nn},'linewidth',2);         %
end         % ------------------------------------------------------------- %
hold off; box on; axis tight;
axis([xlim ylim]);                                                          % rescale the figure
J1=sum(isnan(ZZ)); J2=J1==T; pl(find(J2))=[];  
if lg 
    leg1=legend(pl,legenda,'location','best');             % legend
    set(leg1,'Interpreter','latex');
end
set(gca,'Layer','top');
