%%%% Dynamic Factor Model PhD Course Project
%%%% Pablo BARRIO & Kilian ANDRU
%%%% Estimate Output Gap with Non-Stationnary data

clear
close all
clc

%giorno='20190329';                                              % RESTAT SUBMISSION
%gm=datenum(str2num(giorno(1:4)),...                     % -----------------
%str2num(giorno(5:6)),str2num(giorno(7:8)));         % -----------------
%Filename=['USDB_Haver_' num2str(gm)];                   % file name    
tresh=10^(-2);                                                  % tolerance treshold for EM algorithm
star=10^(-5);                                                   % Initial variance R
maxiter=50;                                                     % max number iteratioin EM algorithm
trans=3;                                                        % data transformation   
out=2;                                                          % outlier treatment
q=4;                                                            % number of factors
s=0;                                                            % lags in the factor loadings
d=q-1;                                                  % number of common cycles
p=3;                                                            % lags VAR
det=1;                                                          % parameters DFM
GDO=0;                                                          % impose GDO restrictions
model='VAR';                                                    % determines law of motion for the factors
m=0;                                                            % parameters Robinson-Yao-Zhang
cc=10;                                                          % obs to exclude because of initial condition
TR1=[];                                      % Variables for which I overwrite the trend test
I0=[];                                  % restrictions for EM algorithm
I1=[];                                                          % -----------------------------
rr=ones(q,1);                                           % -----------------------------
TV.id={1,33};                                                 % time varying parameters
TV.Type={['trend';'none '],'mean'};                             % -----------------------
TV.q0=[10^(-3), 10^(-2)];                                       % initial variance for TV states
nboot=5;                                                     % Number of bootstrap
iter = 100;
method_data = 4;                                                %Method to impute covid data
country = "DE";                                                  %Country
trans_treatment = 'light';                                      %Transformation for data treatment
Block = 9;                                                       %Average size of blocks for bootstrap procedure

[Y, Names, dates] = Data_treatment(country, trans_treatment, method_data );

%%% =============================== %%%
%%%  Initialize the model with PCA  %%%
%%% =============================== %%%
y=ML_diff(Y);                                                           % Data in 1st Differences
[T,N]=size(y);                                                          % size of the panel
[yy, my, sy]=ML_Standardize(y);                                         % Standardize
[f,lambda]=ML_efactors2(yy,q,2);                                        % estimate factor loadings, aka DFM on \Delta y_t
TT=(1:T+1)';                                                            % time trend
X=NaN(T+1,N); bt=X; b=zeros(N,2);                                       % preallocates variables for detrending
J=BL_TestLinearTrend(y); J(TR1)=1;                                      % Identify variables to be detrended
[X(:,J),bt(:,J),b(J,:)]=ML_detrend(Y(:,J));                             % Detrend variables to be detrended
X(:,~J)=Y(:,~J)-repmat(mean(Y(:,~J)),T+1,1);                            % Demean variables not to be detrended
bt(:,~J)=repmat(mean(Y(:,~J)),T+1,1);                                   % ---------------------
b(~J,1)=mean(Y(:,~J))';b(~J,2)=0;                                       % ---------------------
if GDO==1 % ----------------------------------------------------------- % Restrictions for GDO
    lambda(1:2,:)=repmat(mean(lambda(1:2,:)),2,1);                      % same loadings
    b(1:2,2)=mean(b(1:2,2));                                            % same slope, different constant
    bt(:,1:2)=repmat(b(1:2,1)',T+1,1)+TT*b(1:2,2)';                     % same linear trend
    X(:,1:2)=(Y(:,1:2)-bt(:,1:2));                                      % Detrended GDP GDI
    sy(1:2)=repmat(mean(sy(1:2)),1,2);                                  % same std for standardization
end        % -------------------------------------------------------------- %
Z=X./repmat(sy,size(X,1),1);                                            % BLL - Standardization
F=Z*lambda/N;                                                           % Factors in levels as in BLL
[A0,v,AL]=ML_VAR(F,p,1);                                                % Estimate VAR
xi=Z-F*lambda';                                                         % idiosyncratic component
    
%%% ============================ %%%
%%%  Estimate the model with EM  %%%
%%% ============================ %%%
Z2=(Y-repmat(b(:,1)',size(Y,1),1))./repmat(sy,size(Y,1),1);             % BL standardizations: (Y_t - a), i.e. approx. centered levels, divided by sy
b1=b./repmat(sy',1,2); b1(:,1)=0;                                       % divide slope by sy, constant=0;
NSDFM_SS=ML_NSDFM_SS_GDO_TV(AL,lambda,v,q,s,xi,F,rr,I0,I1,b1,A0,TV,star);   % State-space representation
EM=ML_NSDFM_EM_GDO_TV(Z2,NSDFM_SS,iter,tresh,cc,GDO);                   % EM-Algorithm
T2=size(Z2,1); start=1;
BT=TT*EM.b2';                                                           % ML estimates of linear trend
isTV=find(NSDFM_SS.type2==10)+max(p*q,q*(s+1));                         % identifies TV coefficients
for ii=1:length(isTV) % --------------------------------------------------- % TV slopes or means
    id=find(EM.Lambda(:,isTV(ii)));
    BT(:,id)=EM.xitT(:,isTV(ii))*ones(1,length(id));
end                   % --------------------------------------------------- %
FF=EM.xitT(:,1:q); FF1=ML_lag2(FF,s,0);                                 % Store ML Factors
L=EM.Lambda(:,1:(s+1)*q);                                               % Store ML loadings
L2=NaN(N,q,s+1); for ss=1:s+1; L2(:,:,ss)=L(:,(ss-1)*q+1:ss*q); end     % -----------------

%%% =================== %%%
%%%  Common Components  %%%
%%% =================== %%%
T3=length(FF1(cc:end,:));
SY=repmat(sy,T3,1);
MY=repmat(b(:,1)',T3,1);
start2=start+s+cc-1;
Y2=Y(start2:end,:);
chi=(FF1(cc:end,:)*L'+BT(start2:end,:)).*SY+MY;                            % common component
zeta=Y2-chi;  

%%% ================================= %%%
%%%  Estimate common trend and cycles of Factors with EM  %%%
%%% ================================= %%%

disp("Decomposition");

%%% Initialization of EM through PCA %%%

[t,n]=size(FF);                                                          % size of the panel of factors
q2 = 1;                                                                  % number of trends
p2 = 1;                                                                  % number of lags in VAR of trends
[tau00_1 ,Psi0]=ML_efactors2(FF, q2 ,2);                                      % estimate psi aka the loadings of the trends
Z=X./repmat(sy,size(X,1),1);                                             % BLL - Standardization
tau00 = FF(1,:)*Psi0/N;                                                            % trends
A = 1;
%[A0,v,AL]=ML_VAR(tau,p2,1);                                              % Estimate VAR
w00 = FF-tau00*Psi0';                                                         % idiosyncratic component/cycles
R = cov(w00);
Q = var(ML_diff(tau00_1));
P00 = var(tau00_1);

EM2 = EM_decomposition(FF,tau00,P00,Psi0,R,Q,q2,maxiter,tresh,cc);

tautau = EM2.xitT;
Qtau = EM2.Q;
Rww = EM2.R;
Psi = EM2.Psi;
psi_tau = tautau*Psi';                                                  %Common trend
ww = FF - psi_tau;                                                      % cycles


chit = (psi_tau(cc:end,:)*L2'+BT(start2:end,:));                 % common stochastic trend plus linear trend
chist = (psi_tau(cc:end,:)*L2');                                 % common stochastic trend
chic = (ww(cc:end,:)*L2');                                          % common temporary
chint = (psi_tau(cc:end,:)*L2');                                    % common stochastic trend
chilt = (BT(start2:end,:));                                         % linear trend
    

%%% ================================= %%%
%%%  Bootstrap procedure  %%%
%%% ================================= %%%

disp("Bootstrap");

eta=Z2-EM.xitT*EM.Lambda';
chiB=NaN*repmat(chi,1,1,nboot); chicB=chiB; chitB=chiB;
tic
type2=NSDFM_SS.type2;

parfor bb=1:nboot; disp(bb)

    [Z2s,xitTs]=ML_NSDFM_DGP_GDO_TV(T+1,N,q,s,p,EM,eta,type2,[1 0],Block);
    Boot=Bootstrap_GDO_TV(Z2s,NSDFM_SS,GDO,iter,tresh,cc,EM,xitTs,maxiter);

    chiB(:,:,bb)=Boot.chi;
    chicB(:,:,bb)=Boot.chic;
    chitB(:,:,bb)=Boot.chit;

end
toc


 %%% ====================== %%%
 %%% == ---------------- == %%%
 %%% == Plotting Results == %%%
 %%% == ---------------- == %%%
 %%% ====================== %%%

 alpha1=16; alpha2=32;
 Dates = datenum(dates);
 Dates2=Dates(start2:end); Dates2d=Dates2(2:end); Dates4q=Dates2(5:end);     % New Dates for plots
 j0=1;                                                                       % Strarting point for all graphs
 colore={[0 0.45 0.74],[0.64 0.08 0.18],[0.85 0.33 0.1],...                  % colors for gaps blue, red, orange
        [0.93 0.69 0.13],[0 0 0],[0.47 0.67 0.19],[0.49 0.18 0.56]};         % yellow,black, green, purple
 nm=size(chic,2);                                                            % number of possible decompositions
 qsp=[num2str(q) num2str(s) num2str(p)];
 DD=Dates2(j0:end);


LS={'k-','k--','k:','k-.'};
nomefile='Figure_';


lg1={'Estimation'};
lg2={'GDP','Estimation'};

%%% ============ %%%
%%%  Output gap  %%%
%%% ============ %%%

% Output Gap - levels
ZZ=chic(:,1); ZZ=ZZ(j0:end,:);                          
ZZb{1}=BL_Band(chic(j0:end,1),squeeze(chicB(j0:end,1,:)),alpha1);           % -------------------
ZZb{2}=BL_Band(chic(j0:end,1),squeeze(chicB(j0:end,1,:)),alpha2);           % -------------------
SizeXY=ML_SizeXY(DD,ML_minmax([ZZ ZZb{1}]),1,5,[-60 60]);                   % -------------------
BL_ShadowPlotBW(ZZb,ZZ,DD,lg1,SizeXY,LS([1 2]),2,4)
print('-dpng','-vector','-r600',[nomefile 'OG']);     % -------------------
title('Output Gap - Level','fontweight','bold','fontsize',16)   % -------------------
%print('-depsc','-vector','-r600',[nomefile 'OG']);


 % Output Gap - 4Q
ZZ= chic(:,1); ZZ=ML_diff(ZZ(j0:end,:),1);               
ZZb{1}=BL_Band(ZZ,squeeze(ML_diff(chicB(j0:end,1,:),1)),alpha1);       % ---------------
ZZb{2}=BL_Band(ZZ,squeeze(ML_diff(chicB(j0:end,1,:),1)),alpha2);       % ---------------
SizeXY=ML_SizeXY(Dates4q(j0:end),ML_minmax([ZZ ZZb{1}]),1,5,[-60 60]);      % ---------------
BL_ShadowPlotBW(ZZb,ZZ,Dates4q(j0:end),lg1,SizeXY,LS([1 2]),2,4)
print('-dpng','-vector','-r600',[nomefile 'OG4Q']);    % ---------------
title('Output Gap - 4Q % changes','fontweight','bold','fontsize',16) 


chit2 = chit.*SY+MY;
chitB2 = chitB.*SY+MY;
chit3 = exp(chit2)./1000;
chitB3 = exp(chitB2)./1000;

Y3 = exp(Y2)./1000;

% Potential output (levels)
YY=[Y3(:,1)  chit3(:,1)]; YY=YY(j0:end,:);                       
YYb{1}=BL_Band(YY(:,2),squeeze(chitB3(j0:end,1,:)),alpha1);                         % ------------------------
YYb{2}=BL_Band(YY(:,2),squeeze(chitB3(j0:end,1,:)),alpha2);                         % ------------------------
SizeXY=ML_SizeXY(DD(1:60),ML_minmax([YY YYb{1}]),1,5);
YYY{1}=YYb{1}(1:60,:); YYY{2}=YYb{2}(1:60,:);                                   % ------------------------
BL_ShadowPlotBW(YYY,YY(1:60,:),DD(1:60),lg2,SizeXY,LS([2 1]),2,4,50)
print('-dpng','-vector','-r600',[nomefile 'PO_7084']);	% ------------------------
title('Potential Output - 4Q % changes','fontweight','bold','fontsize',16) 


% Potential Output - 4Q
ZZ=[Y3(:,1) chit3(:,1)]; ZZ=ML_diff(ZZ(j0:end,:),1);        
ZZb{1}=BL_Band(ZZ(:,2),squeeze(ML_diff(chitB3(j0:end,1,:),1)),alpha1);       % ---------------------
ZZb{2}=BL_Band(ZZ(:,2),squeeze(ML_diff(chitB3(j0:end,1,:),1)),alpha2);       % ---------------------
SizeXY=ML_SizeXY(Dates4q(j0:end),ML_minmax([ZZ ZZb{1}]),1,5,[-60 60]);      % ---------------------
BL_ShadowPlotBW(ZZb,ZZ,Dates4q(j0:end),lg2,SizeXY,LS([2 1]),2,4,20)
print('-dpng','-vector','-r600',[nomefile 'PO']);      % ---------------------
title('Potential Output - 4Q % changes','fontweight','bold','fontsize',16) 	% ---------------------


% Contribution to GDP growth (BL)
 ZZ=ML_diff([chit(:,1) chic(:,1) zeta(:,1)],1); YY=ML_diff(Y2(:,1),1).*100;               
ML_ContributionGraphBW(YY(j0:4:end),ZZ(j0:4:end,:),Dates4q(j0:4:end),...     %
    {'GDP','Potential output','Output gap','Idiosyncratic'});                       % --------------------------
print('-dpng','-vector','-r600',[nomefile 'DecompGDP_BL']);   % --------------------------
title('GDP growth decomposition','fontweight','bold','fontsize',16)

Dates_plot = year(Dates)+ (month(Dates)-1)/11;

% Trend - levels
ZZ=[tautau zeros(length(Dates),1)]; ZZ=ZZ(j0:end,:); % -------------------
plot(Dates_plot(j0:end,:),ZZ);
print('-dpng','-vector','-r600',[nomefile 'trend']);     % -------------------
title('Trend - Level','fontweight','bold','fontsize',16)   % -------------------


%  % Trend - 4Q
ZZ= [diff(tautau(:,1)) zeros(length(Dates)-1,1)];ZZ=ZZ(j0:end,:);  % -------------------
plot(Dates_plot(j0+1:end),ZZ);
print('-dpng','-vector','-r600',[nomefile 'trend4q']);     % -------------------
title('Trend - First difference','fontweight','bold','fontsize',16)   % -------------------


for plot_ii = 1:q
%  cycle - levels
ZZ=[ww(:,plot_ii) zeros(length(Dates),1)]; ZZ=ZZ(j0:end,:); % -------------------
plot(Dates_plot(j0:end,:),ZZ);
print('-dpng','-vector','-r600',[nomefile 'cycle_' num2str(plot_ii)]);     % -------------------
title(['Cycle ' num2str(plot_ii) ' - Level'],'fontweight','bold','fontsize',16)   % -------------------


%  cycle - 4Q
ZZ= [diff(ww(:,plot_ii)) zeros(length(Dates)-1,1)];ZZ=ZZ(j0:end,:);  % -------------------
plot(Dates_plot(j0+1:end),ZZ);
print('-dpng','-vector','-r600',[nomefile 'cycle4q'  num2str(plot_ii)]);     % -------------------
title(['Cycle ' num2str(plot_ii) ' - Fist diffrence'],'fontweight','bold','fontsize',16)   % -------------------
end


