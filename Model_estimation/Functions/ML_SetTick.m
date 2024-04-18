

function [TICK,xtl,J]=ML_SetTick(Dates,ny,mq)

temp=datevec(Dates); TICK=find(temp(:,2)==mq);                               % Tick at each mq month/quarter of the year
TICK2=TICK(mod(temp(TICK,1),ny)==0);                                         % xtick every ny years
xtl=cellstr(repmat(' ',length(TICK),1));                                    % define the stile of the ticklabel
J=mod(temp(TICK,1),ny)==0;
xtl(J)=cellstr(datestr(Dates(TICK2),'yyyy'));          % ---------------------------------
