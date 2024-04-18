% BL_PlotRealTime2 - (ZZ,Dates,Titolo,legenda,asseX,sizeXY,ng,colore)
function BL_PlotRealTimeBW(ZZ,Dates,Titolo,legenda,asseX,sizeXY,ng,ny,mq,ypace)

if nargin<7; ng=0; end
if nargin<6; sizeXY=[]; end
[T,N]=size(ZZ);    


try isnan(ny); catch, ny=5; end
try isnan(mq); catch, mq=3; end
try isnan(ypace);
    ndec=['%.' num2str(ML_ndecimal(ypace)) 'f'];
catch, ypace=1; ndec='%.0f';
end
grey=[.8 .8 .8];                                                            % grey is for grid
lavender=[230,230,250]/255;                                                 % lavender is for shaded area


if isempty(sizeXY);                                                         % if size of figure not provided
    xlim=Dates([1 end])'; ylim=[ML_min(ZZ) ML_min(ZZ,2)];                   % ------------------------------
else xlim=sizeXY(1:2); ylim=sizeXY(3:4);                                    % size of figure provided by User
end                                                                         % ------------------------------

temp=datevec(Dates); TICK=find(temp(:,2)==mq);                               % Tick at each first quarter of the year
TICK2=TICK(mod(temp(TICK,1),ny)==0);                                         % xtick every five years
xtl=cellstr(repmat(' ',length(TICK),1));                                        % define the stile of the ticklabel
xtl(mod(temp(TICK,1),ny)==0)=cellstr(datestr(Dates(TICK2),'yyyy'));          % ---------------------------------
if sum(sign(ylim))==0;
    t1=0:-ypace:ML_round(ylim(1),0); t2=0:ypace:ML_round(ylim(2),0);
    ytick=sort([t1(2:end) t2]); 
elseif sum(sign(ylim))>0;
    ytick=ML_round(ylim(1),0):ypace:ML_round(ylim(2),0);       
else yytick=ML_round(ylim(1),0):ypace:ML_round(ylim(2),0);
end

load NBER_Recessions; NBER2=NBER(:,2)*ylim;


        %%% ======================= %%%
        %%%         Graphing        %%%
        %%% ======================= %%%
        
axes('Parent',figure,'FontSize',10); ML_FigureSize,hold on;
ha=area(NBER(:,1),[NBER2(:,1) NBER2(:,2)-NBER2(:,1)],'linestyle','none');       % NBER Recessions
set(ha(1), 'FaceColor', 'none'); set(ha(2), 'FaceColor', [.9 .9 .9])              % ---------------
axis([xlim ylim])                                                               % define size of the figure
% K=cat(1,find(ML_diff(NBER(:,2))==1)+1,find(ML_diff(NBER(:,2))==-1));
% for ii=1:length(K); line(repmat(NBER(K(ii),1),1,2), NBER2(K(ii),:),'color','k'); end
set(gca,'Xtick', Dates(TICK),'Xticklabel',xtl);                                 % set xtick
set(gca,'Ytick',ytick,'Yticklabel',num2str(ytick',ndec));
plot(Dates,zeros(T,1),'k','linewidth',.5);                                      % line at zero

plot(Dates,ZZ(:,2:N-2),'color',[.7 .7 .7],'linewidth',.5);         
pl(1)=plot(Dates,ZZ(:,1),'k--','linewidth',2);      
pl(2)=plot(Dates,ZZ(:,N-1),'k','linewidth',2);      
pl(3)=plot(Dates,ZZ(:,N),'k.','markersize',15);      

hold off; box on; axis tight;
axis([xlim ylim]);                                                          % rescale the figure  
if isempty(legenda); else                                                   % legend
    J1=sum(isnan(ZZ)); J2=J1==T; pl(find(J2))=[];
    legend(pl,legenda,'location','best'); 
end             
if isempty(Titolo); else title(Titolo,'fontsize',18,'fontweight','bold'); end   % title
if isempty(asseX); else xlabel(asseX,'fontsize',14,'fontweight','bold'); end    % xlabel
% ML_EliminateUpperRightTicks;                                                    % eliminate upper ticks
set(gca,'Layer','top');
