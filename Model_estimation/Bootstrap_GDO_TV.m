function Boot=Bootstrap_GDO_TV(Z2,NSDFM_SS,GDO,iter,tresh,cc,EM,xitTs,maxiter)

s=NSDFM_SS.s; p=NSDFM_SS.p; q=NSDFM_SS.q;

xitT=EM.xitT;

NSDFM_SS.P00=EM.PtT(:,:,1);
NSDFM_SS.F00=EM.xitT(1,:)';


%%% ================================= %%%
%%%  Durbin Koopman 2002 Algorithm 2  %%%
%%% ================================= %%%
type2=NSDFM_SS.type2;
EM=ML_NSDFM_EM_GDO_TV_sim(Z2,NSDFM_SS,iter,tresh,cc,GDO);                       % EM-Algorithm                      
xitTb=xitT-EM.xitT+xitTs;                                                   % Centering the states

[T2, N]=size(Z2); start=1;
TT=(1:T2)';                                                                 % time trend
BT=TT*EM.b2';                                                               % ML estimates of linera trend
isTV=find(type2==10)+max(p*q,q*(s+1));                                      % identifies TV coefficients
for ii=1:length(isTV) % --------------------------------------------------- % TV slopes or means
    id=find(EM.Lambda(:,isTV(ii)));
    BT(:,id)=xitTb(:,isTV(ii))*ones(1,length(id));                     
end                   % --------------------------------------------------- %
FF=xitTb(:,1:q); FF1=[]; for ss=0:s;FF1=cat(2,FF1,FF(s-ss+1:end-ss,:));end   % Store ML Factors
for ss=1:s+1; lambda2(:,:,ss)=EM.Lambda(:,(ss-1)*q+1:ss*q); end                % Store ML loadings
L=reshape(lambda2,N,q*(s+1));
L2=NaN(N,q,s+1); for ss=1:s+1; L2(:,:,ss)=L(:,(ss-1)*q+1:ss*q); end

%%% ================= %%%
%%%  Estimate trend   %%%
%%% ================= %%%
%%% Initialization of EM through PCA %%%

q2 = 1;                                                                  % number of trends
[tau00_1 ,Psi0]=ML_efactors2(FF, q2 ,2);                                      % estimate psi aka the loadings of the trends
tau00 = FF(1,:)*Psi0/N;                                                            % trends
%[A0,v,AL]=ML_VAR(tau,p2,1);                                              % Estimate VAR
w00 = FF-tau00*Psi0';                                                         % idiosyncratic component/cycles
R = cov(w00);
Q = var(ML_diff(tau00_1));
P00 = var(tau00_1);

EM2 = EM_decomposition(FF,tau00,P00,Psi0,R,Q,q2,maxiter,tresh,cc);

tautau = EM2.xitT;
Psi = EM2.Psi;
psi_tau = tautau*Psi';                                                  %Common trend
ww = FF - psi_tau;                                                      % cycles

start2=start+s+cc-1;
chi=FF1(cc:end,:)*L'+BT(start2:end,:);                                  % common component
chit = (psi_tau(cc:end,:)*L2'+BT(start2:end,:));                        % common stochastic trend plus linear trend                               
chic = (ww(cc:end,:)*L2');                                              % common temporary
chilt = (BT(start2:end,:));                                             % linear trend


eta=Z2(start2:end,:)-xitTb(cc:end,:)*EM.Lambda';                            % Stationary residuals (aka xi~I(0) plus near zero)
zeta=0*eta;                                                                 % xi~I(1)
isxi=find(type2==1)+12;                                                     % identifies states that are idiosyncratic components
istrend=find(EM.mu~=0); istrend(1:q)=[];                                       % identifies states with linear trend
isc2=isxi(~ismember(isxi,istrend));                                         % identifies states that are idio RWnLT
isc4=isxi(ismember(isxi,istrend));                                          % identifies states that are idio RWwLT
J2=find(sum(EM.Lambda(:,isc2),2)==1); nJ2=ismember(1:N,J2);                    % identifies variables that have idio RWnLT 
J4=find(sum(EM.Lambda(:,isc4),2)==1); nJ4=ismember(1:N,J4);                    % identifies variables that have idio RWwLT 
zeta(:,nJ2)=xitTb(cc:end,isc2);                                           % idio for RWnLT
zeta(:,nJ4)=xitTb(cc:end,isc4)-chilt(:,nJ4);                              % idio for RWwLT

Boot.chi=chi;
Boot.chic=chic;
Boot.chit=chit;
Boot.zeta=zeta;
Boot.eta=eta;
Boot.chilt=chilt;
Boot.xitTb=xitTb;
Boot.LambdaS=EM.Lambda;