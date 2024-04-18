
% Algorithm 2 Durbin Koopman 2002
% xitTb = \tilde{\alpha}
% xitT = \hat{\alpha}
% EM.xitT = \hat{\alpha}^\dag
% xitTs = \alpha^\dag

function Boot=ML_NSDFM_Boot_GDO_TV(Z2,NSDFM_SS,GDO,iter,tresh,cc,EM,xitTs,d)

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

%%% ================================= %%%
%%%  Estimate trend with PCA on LRCV  %%%
%%% ================================= %%%
[V, D]= eig(cov(FF(cc:end,:)));                                             % eigenvalue and eigenvenctors of VCV
[~, t2]=sort(diag(D),'descend'); V=V(:,t2); D=D(t2,t2);                     % sort eigenvalues and eigenvectors
FFt = FF(cc:end,:)*V(:,1:q-d);                                              % TREND Factorss
FFc = FF(cc:end,:)*V(:,q-d+1:q);                                            % CYCLE Factors
FFt1=[]; FFc1=[];   
for ss=1:s+1;     
    FFt1=cat(2,FFt1,FFt(s-ss+2:end-ss+1,:));                                % lagged trend factors
    FFc1=cat(2,FFc1,FFc(s-ss+2:end-ss+1,:));                                % lagged cycle factors
    lambdat(:,(ss-1)*(q-d)+1:ss*(q-d)) = lambda2(:,:,ss)*V(:,1:q-d);        % trend loadings
    lambdac(:,(ss-1)*d+1:ss*d) = lambda2(:,:,ss)*V(:,q-d+1:q);              % cycle loadings    
end

start2=start+s+cc-1;
chi=FF1(cc:end,:)*L'+BT(start2:end,:);                                    % common component
chit = (FFt1*lambdat'+BT(start2:end,:));                                    % common stochastic trend plus linear trend
chic = (FFc1*lambdac');                                                     % common temporary
chilt=BT(start2:end,:);

eta=Z2(start2:end,:)-xitTb(cc+1:end,:)*EM.Lambda';                            % Stationary residuals (aka xi~I(0) plus near zero)
zeta=0*eta;                                                                 % xi~I(1)
isxi=find(type2==1)+12;                                                     % identifies states that are idiosyncratic components
istrend=find(EM.mu~=0); istrend(1:q)=[];                                       % identifies states with linear trend
isc2=isxi(~ismember(isxi,istrend));                                         % identifies states that are idio RWnLT
isc4=isxi(ismember(isxi,istrend));                                          % identifies states that are idio RWwLT
J2=find(sum(EM.Lambda(:,isc2),2)==1); nJ2=ismember(1:N,J2);                    % identifies variables that have idio RWnLT 
J4=find(sum(EM.Lambda(:,isc4),2)==1); nJ4=ismember(1:N,J4);                    % identifies variables that have idio RWwLT 
zeta(:,nJ2)=xitTb(cc+1:end,isc2);                                           % idio for RWnLT
zeta(:,nJ4)=xitTb(cc+1:end,isc4)-chilt(:,nJ4);                              % idio for RWwLT


Boot.chi=chi;
Boot.chic=chic;
Boot.chit=chit;
Boot.zeta=zeta;
Boot.eta=eta;
Boot.chilt=chilt;
Boot.xitTb=xitTb;
Boot.LambdaS=EM.Lambda;


%%% Rodriguez Ruiz CSDA 2012. TO BE CLEANED
% [xitTd,~,~,~,~,~,~,As,LambdaS,Rs,Qs,b2,mus]=ML_DFM_EM_GDO_TV_sim...           % EM-Algorithm
%     (Z2,F00,P0s,A,Lambda,R,Q,s,q,d,p,mu,type2,iter,tresh,cc,GDO);  % ------------
% [xitt,xittm,Ptt,Pttm]=ML_KalmanFilter2(F00,P0s,Z2,As,LambdaS,Rs,Qs,mus);   % Kalman Filter                           
% [xitTb]=ML_KalmanSmoother2(A,xitt',xittm',Ptt,Pttm,LambdaS,Rs);       % Kalman Smoother