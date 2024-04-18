
function BL_ShadowPlotBW(XX,ZZ,Dates,legenda,sizeXY,LS,ny,mq,ypace,ytickx,LW)



[T,N]=size(ZZ);    


grey={[.8 .8 .8], [.7 .7 .7]};                                                            % color for grid (grey)
cww=[.9 .9 .9];                                                      % color for NBER recession (lavender)

try isnumeric(ny); if isempty(ny); ny=5; end;  catch; ny=5; end
try isnan(mq); if isempty(mq); mq=3; end; catch; mq=3; end
try isempty(ypace); ndec=['%.' num2str(ML_ndecimal(ypace)) 'f'];
    if isempty(ypace); ypace=1; ndec='%.0f'; end
catch, ypace=1; ndec='%.0f';
end

try isnan(LW); if isempty(LW); LW=2; end; catch; LW=2; end

if isempty(sizeXY)                                                          % if size of figure not provided
    xlim=Dates([1 end])'; ylim=[ML_min(ZZ) ML_min(ZZ,2)];                   % ------------------------------
else; xlim=sizeXY(1:2); ylim=sizeXY(3:4);                                    % size of figure provided by User
end                                                                         % ------------------------------

if sum(sign(ylim))==0  % -------------------------------------------------- % Set-up yticks
    t1=0:-ypace:ML_round(ylim(1),0); t2=0:ypace:ML_round(ylim(2),0);
    ytick=sort([t1(2:end) t2]); 
elseif sum(sign(ylim))>0
    ytick=ML_round(ylim(1),0):ypace:ML_round(ylim(2),0);       
else; ytick=ML_round(ylim(1),0):ypace:ML_round(ylim(2),0);
end

ytick(ismember(ytick,ylim))=[];


% [peak, trough]=ML_NBERrecession;                                            % upload data for NBER Recessions
% for jj=1:length(peak); NBER{jj}=[peak(jj) trough(jj)]; end                  % build the variable for the graph

temp=datevec(Dates); TICK=find(temp(:,2)==mq);                              % Tick at each mq of the year
TICK2=TICK(mod(temp(TICK,1),ny)==0);                                        % grid every ny years
xtl=cellstr(repmat(' ',length(TICK),1));                                    % define the stile of the ticklabel
xtl(mod(temp(TICK,1),ny)==0)=cellstr(datestr(Dates(TICK2),'yyyy'));         % ---------------------------------

load NBER_Recessions; NBER2=NBER(:,2)*ylim;

if iscell(XX); J=length(XX); WW=XX; else; J=1; WW{1}=XX; end

        %%% ======================= %%%
        %%%         Graphing        %%%
        %%% ======================= %%%
        
axes('Parent',figure,'FontSize',10); ML_FigureSize,hold on;
ha=area(NBER(:,1),[NBER2(:,1) NBER2(:,2)-NBER2(:,1)],'linestyle','none');   % NBER Recessions
set(ha(1), 'FaceColor', 'none'); set(ha(2), 'FaceColor', cww)               % ---------------
axis([xlim ylim])                                                           % define size of the figure
set(gca,'Xtick', Dates(TICK),'Xticklabel',xtl);                             % set xtick

try isnumeric(ytickx);
    if isnumeric(ytickx)
        yX=repmat(' ',length(ytick),ytickx);                                % add blank spaces to y-tick
        set(gca,'Ytick',ytick,'Yticklabel',[num2str(ytick',ndec) yX]);      % set ytick
    else
        set(gca,'Ytick',ytick,'Yticklabel',num2str(ytick',ndec));           % set ytick
    end
catch
    set(gca,'Ytick',ytick,'Yticklabel',num2str(ytick',ndec));               % set ytick
end

for jj=1:J
    ha=area(Dates,[WW{jj}(:,2) WW{jj}(:,1)-WW{jj}(:,2)],'linestyle','none');
    set(ha(1), 'FaceColor', 'none'); set(ha(2), 'FaceColor', grey{jj})
end
plot(Dates,zeros(T,1),'k','linewidth',.5);                       % line at zero

for nn=1:N  % ------------------------------------------------------------- % plot each variable separately
    pl(nn)=plot(Dates,ZZ(:,nn),LS{nn},'linewidth',LW);          % 
end         % ------------------------------------------------------------- %
hold off; box on; axis tight;
axis([xlim ylim]);                                                              % rescale the figure
J1=sum(isnan(ZZ)); J2=J1==T; pl(find(J2))=[];  
if isempty(legenda); else legend(pl,legenda,'location','best','interpreter','latex'); end             % legend
set(gca,'Layer','top');