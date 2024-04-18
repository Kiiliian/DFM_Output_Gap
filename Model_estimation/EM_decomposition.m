 % [xitT,PtT,PtTm,xitt,xittm,Ptt,Pttm,A,Psi,R,Q,b1,VECM]=...
%       EM_decomposition(Y,F00,P0s,A,Psi,R,Q,s,q,d,p,mu,type,model,maxiter,tresh,cc)
% 
%       Y - data with deterministic component
%     F00 - initial values for the states
%     P0s - initial variance for the states
%       A - VAR(1) for the states
%  Psi - Loadings
%       R - covariance matrix errors
%       Q - covariance matrix shocks
%       s - number of lags in the factor loadings
%       q - number of shocks
%       p - number of lags VAR factors
%      mu - constant in Factors VAR, and slope for LT
%    type - identify type of state (TV,RW,WNwLT,RWwLT)
%   model - {'VAR','VECM'}
% maxiter - max number iterations
%   tresh - treshold for stopping algorithm
%      cc - # of obs eliminated before estimateing parameters
% 

function EM=EM_decomposition(Y,tau00,P00,Psi,R,Q,q,maxiter,tresh,cc)
% function EM=EM_decomposition(Y,F00,P0s,A,Psi,R,Q,s,q,d,p,mu,type,model,maxiter,tresh,cc)

A=1;
s=0;
p=1;


[T, N]=size(Y);                                                             % size of the panel
s1=s+1;                                                                     % number of lags plus one 
r=q*(s1);                                                                   % number of "static" factors
pq=q*p;                                                                     % Number of states in the VAR                                                                 % time trend                                                  % TV deterministic components

%%% =================== %%%
%%%  Start with E-Step  %%%
%%% =================== %%%
[xitt,xittm,Ptt,Pttm,loglik]=ML_KalmanFilter2(tau00,P00,Y,A,Psi,R,Q);   % Kalman Filter                           
[xitT,PtT,PtTm]=ML_KalmanSmoother2(A,xitt',xittm',Ptt,Pttm,Psi,R);       % Kalman Smoother
tau00=xitT(1,:)'; P00=PtT(:,:,1);                                             % initial conditions
tau=xitT(:,1:q);                                                              % factors
tautau=[];

for ss=0:s
    tautau=cat(2,tautau,tau(s-ss+1:end-ss,:));
end                                                                        % lagged factors 0:s

loglik1=loglik;                                                             % likelihood

disp('Running EM algorithm: ')

disp(['0: ' num2str(loglik)])
for jj=1:maxiter         
    
    %%%% ======================== %%%%
    %%%% ====                ==== %%%%
    %%%% ====    M - Step    ==== %%%%
    %%%% ====                ==== %%%%
    %%%% ======================== %%%%
    
    %%% =========================== %%%
    %%% 1: Estimate factor Loadings %%%
    %%% =========================== %%%
    yy=ML_center(Y)';                                                       % endogenous variables,     
    xx=cat(2,zeros(r,s),ML_center(tautau)');                                    % exogenous variables, aka the factors    
    num=zeros(N,r); den=zeros(r,r);                                         % preallocates
    for tt=cc+1:T  % ------------------------------------------------------ % \sum_{t=1}^T
        num=num+(yy(:,tt))*xx(:,tt)';                                       % (y_t - \xi_t) F_{t|T}'
        den=den+xx(:,tt)*xx(:,tt)'+PtT(1:r,1:r,tt);                         % F_{t|T}F_{t|T}'+P_{t|T}
    end    
    L=num/den;                                                              % factor loadings       
    Psi(1:N,1:r) = L;                                                    % store factor loadings


    %%% ======================================================= %%%
    %%% 2: Estimate parameters for law of motion of the Factors %%%
    %%% ======================================================= %%%    
    yy=cat(2,zeros(q,1),tau(2:T,:)');                                         % F_t
    xx=cat(1,ones(1,T),cat(2,ML_lag(xitT(:,1:pq),1)',zeros(pq,1)));         % F_{t-1}
    EF1=zeros(q,pq+1); EF1F1=zeros(pq+1,pq+1); EF=zeros(q,q);               % initialize
    for tt=cc+2:T % ------------------------------------------------------ % \sum_{t=1}^T
        EF1=EF1+yy(:,tt)*xx(:,tt-1)'+[zeros(q,1) PtTm(1:q,1:pq,tt)];        % E(F_t F_{t-1})
        EF1F1=EF1F1+xx(:,tt-1)*xx(:,tt-1)'+blkdiag(0,PtT(1:pq,1:pq,tt-1));  % E(F_{t-1} F_{t-1})
        EF=EF+yy(:,tt)*yy(:,tt)'+PtT(1:q,1:q,tt);                           % E(F_t F_t)
    end
    temp=EF1/EF1F1;
    Q(1:q,1:q) = (EF - temp*EF1') / (T-cc);                                 % covariance
        
              
    
    %%% ============================== %%%
    %%% 3: Covariance Prediction Error %%%
    %%% ============================== %%%
    yy=Y';
    xx=xitT';
    R=zeros(N,N);                                                           % preallocates
    for tt=cc+1:T  % ------------------------------------------------------ % \sum_{t=1}^T
        eta=yy(:,tt)-Psi*xx(:,tt);                                       % prediction error
        R = R+ eta*eta'+ Psi*PtT(:,:,tt)*Psi';                        % E(eta_t eta_t')
    end            % ------------------------------------------------------ %
    R = diag(diag(R/(T-cc)));
    
    %%%% ======================== %%%%
    %%%% ====                ==== %%%%
    %%%% ====    E - Step    ==== %%%%
    %%%% ====                ==== %%%%
    %%%% ======================== %%%%            
    [xitt,xittm,Ptt,Pttm,loglik]=ML_KalmanFilter2(tau00,P00,Y,A,Psi,R,Q); % Kalman Filter
    [xitT,PtT,PtTm]=ML_KalmanSmoother2(A,xitt',xittm',Ptt,Pttm,Psi,R);   % Kalman SMoother           
    tau=xitT(:,1:q);                                                          % jj-th step factors    
    tautau=[]; for ss=0:s; tautau=cat(2,tautau,tau(s-ss+1:end-ss,:)); end                 % jj-th step lagged factors 0:s    
    tau00=xitT(1,:)'; P00=PtT(:,:,1);                                         % initial conditions    

    %%%% ================================= %%%%
    %%%% ====                         ==== %%%%
    %%%% ====    Check convergence    ==== %%%%
    %%%% ====                         ==== %%%%
    %%%% ================================= %%%%
    delta_loglik = loglik - loglik1;                                        % logL(t) - logL(t-1)
    avg_loglik = (abs(loglik) + abs(loglik1) + 10^(-3))/2;                  % avg = (|logL(t)| + |logL(t-1)|)/2
    disp([num2str(jj) ': ' num2str(loglik) '  ' ...
        num2str(delta_loglik / avg_loglik)])
    if (abs(delta_loglik) / avg_loglik) < tresh; break; end                 % convergence if |f(t) - f(t-1)| / avg < threshold
    loglik1=loglik;                                                         % store log-likelihood    
end

% stop
if jj==maxiter 
    disp('not converged'); 
else
    disp(['converged after ' num2str(jj) ' iterations'])
end

EM.xitT=xitT;
EM.PtT=PtT;
EM.PtTm=PtTm;
EM.xitt=xitt;
EM.xittm=xittm;
EM.Ptt=Ptt;
EM.Pttm=Pttm;
EM.A=A;
EM.Psi=Psi;
EM.R=R;
EM.Q=Q;

