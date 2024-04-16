%%%% Dynamic Factor Model PhD Course Project
%%%% Pablo BARRIO & Kilian ANDRU
%%%% Estimate Output Gap with Non-Stationnary data

clear all;
clear clc;

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
GDO=1;                                                          % impose GDO restrictions
model='VAR';                                                    % determines law of motion for the factors
m=0;                                                            % parameters Robinson-Yao-Zhang
cc=10;                                                          % obs to exclude because of initial condition
TR1=[6 15 72 73 76:80 87];                                      % Variables for which I overwrite the trend test
I0=[1:2 36 38 40 43 68:71 75];                                  % restrictions for EM algorithm
I1=[];                                                          % -----------------------------
rr=ones(q,1);                                           % -----------------------------
TV.id={1:2,75};                                                 % time varying parameters
TV.Type={['trend';'none '],'mean'};                             % -----------------------
TV.q0=[10^(-3), 10^(-2)];                                       % initial variance for TV states
nboot=1000;                                                     % Number of bootstrap
iter = 50;

[Y, Names, dates, info] = Data_treatment('FR', 'light', 2 );

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

%%% Initialization of EM through PCA %%%

[t,n]=size(FF);                                                          % size of the panel of factors
q2 = 1;                                                                  % number of trends
p2 = 1;                                                                  % number of lags in VAR of trends
[tau00 ,Psi0]=ML_efactors2(FF, q2 ,2);                                      % estimate psi aka the loadings of the trends
Z=X./repmat(sy,size(X,1),1);                                             % BLL - Standardization
%tau00=FF*Psi/N;                                                            % trends
A = 1;
%[A0,v,AL]=ML_VAR(tau,p2,1);                                              % Estimate VAR
w00 = FF-tau00*Psi0';                                                         % idiosyncratic component/cycles
R = cov(w00);
Q = cov(ML_diff(tau00));
P00 = ;

EM2 = EM_decomposition(FF,tau00,P00,Psi0,R,Q,q2,maxiter,tresh,cc);

tautau = EM2.tau;
ww = EM2. ; %
Psi = EM2.Psi;



