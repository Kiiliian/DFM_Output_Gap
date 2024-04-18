% contribution graph for 4-quarter changes
% Y variable to be plotted
% X contribution variables 

function ML_ContributionGraphBW(Y,X,Dates,legenda)

uno=ones(1,3);
cb=[.5*uno; .7*uno;.9*uno];
cww=[.9 .9 .9];                                                      % color for NBER recession (lavender)
grigio=.4*uno;

% try j0=GP.j0; catch; j0=1; end
if ~isempty(legenda); lg=1; else lg=0; end
XX=cat(3,X,X);                                                              % Matrix of negative and positive contribution
for ii=1:size(X,2); XX(X(:,ii)<0,ii,1)=0; XX(X(:,ii)>0,ii,2)=0; end         % --------------------------------------------

try SizeXY=GP.SizeXY; 
catch
    a=mean(ML_diff(Dates));
    SizeXY=ML_SizeXY(Dates,[Y squeeze(sum(XX,2))],1,5); 
    SizeXY(1:2)=SizeXY(1:2)+[-a a];
end

% load NBER_Recessions;         

axes('Parent',figure,'FontSize',12); ML_FigureSize,hold on;
% area(NBER(:,1),SizeXY(3)*NBER(:,2),'linestyle','none' ,'FaceColor', cww);                      % upper shaded area
% area(NBER(:,1),SizeXY(4)*NBER(:,2),'linestyle','none', 'FaceColor', cww);                      % lower shaded area
[TICK2,xtl,J]=ML_SetTick(Dates(1:end),5,3);
set(gca,'Xtick', Dates(TICK2));
axis(SizeXY); 
% gridxy2(Dates(TICK2(J)),get(gca,'ytick'),'color',[.8 .8 .8],'linewidth',1);
bar1 = bar(Dates,XX(:,:,1),'stacked');
bar2 = bar(Dates,XX(:,:,2),'stacked');
for jj=1:size(X,2)
    set(bar1(jj),'FaceColor',cb(jj,:),'edgecolor',grigio);
    set(bar2(jj),'FaceColor',cb(jj,:),'edgecolor',grigio);    
end
line(get(gca,'xlim'),[0 0],'color','k','linewidth',1);
pl1= plot(Dates,Y,'-sk','linewidth',2,'MarkerFaceColor','k');
hold off; box on; set(gca,'Layer','top');
set(gca,'Xticklabel',xtl,'yticklabel',num2str(get(gca,'ytick')','%.0f'));   % set xtick
if lg; legend([pl1 bar1],legenda,'location','SE'); end             % legend
