% ML_efactorML4c - ML estimation I0 DFM, straight EM no initialization, P00=I (for simulations)
% 
% [xitT,PtT,PtTm,xitt,xittm,Ptt,Pttm,A,Lambda,R,Q]=...
%     ML_efactorML4(y,F0,P0,A,Lambda,R,Q,q,p,maxiter,tresh,cc)% 
% 
%       Y - data with deterministic component
%      F0 - initial values for the states
%      P0 - initial variance for the states
%       A - VAR(1) for the states
%  Lambda - Loadings
%       R - covariance matrix errors
%       Q - covariance matrix shocks
%       r - number of factors
%       q - number of shocks
%       p - number of lags VAR factors
% maxiter - max number iterations
%   tresh - treshold for stopping algorithm
%      cc - # of obs eliminated before estimateing parameters
% 


function LLLM=BL_LLLM(y,SS,maxiter,tresh,cc,LAMBDA,CU,ERRE)

try isempty(LAMBDA);
    if isempty(LAMBDA); sl=1; end
    sl=0; 
catch
    sl=1; 
end

try isempty(CU);
    if isempty(CU); sq=1; end
    sq=0; 
catch
    sq=1; 
end

try isempty(ERRE);
    if isempty(ERRE); sr=1; end
    sr=0; 
catch
    sr=1; 
end

A=SS.A;
Lambda=SS.C;
Q=SS.Q;  
P0=SS.p0;
F0=SS.f0;
R=SS.R;   

r=size(Lambda,2);
OPTS.disp = 0;
[T, N]=size(y);                                                             % size of the panel


%%% =================== %%%
%%%  Start with E-Step  %%%
%%% =================== %%%
[xitt,xittm,Ptt,Pttm,loglik]=ML_KalmanFilter2(F0,P0,y,A,Lambda,R,Q);        % Kalman Filter
[xitT,PtT,PtTm]=ML_KalmanSmoother2(A,xitt',xittm',Ptt,Pttm,Lambda,R);       % Kalman Smoother
F00=xitT(1,:)'; P00=PtT(:,:,1);                                             % initial conditions  
F=xitT;                                                                     % factors
loglik1=loglik;                                                             % likelihood

disp([num2str(0) ': ' num2str(loglik)])

for jj=1:maxiter         
    
    %%%% ======================== %%%%
    %%%% ====                ==== %%%%
    %%%% ====    M - Step    ==== %%%%
    %%%% ====                ==== %%%%
    %%%% ======================== %%%%
    
    %%% =========================== %%%
    %%% 1: Estimate factor Loadings %%%
    %%% =========================== %%%
    if sl==1
        yy=ML_center(y)';                                                       % endogenous variables,
        xx=ML_center(F(:,1))';                                                  % exogenous variables, aka the factors
        num=zeros(N,1); den=0;                                                  % preallocates
        for tt=cc+1:T  % ------------------------------------------------------ % \sum_{t=1}^T
            num=num+yy(:,tt)*xx(:,tt)';                                         % y_t F_{t|T}'
            den=den+xx(:,tt)*xx(:,tt)'+PtT(1,1,tt);                         % F_{t|T}F_{t|T}'+P_{t|T}
        end
        Lambda(1:N,1) = num/den;                                                    % store factor loadings
    end
    
    %%% ==================================== %%%
    %%% 2: Estimate covariance of the States %%%
    %%% ==================================== %%%      
    if sq==1
        yy=cat(2,zeros(r,1),ML_center(F(2:T,:))');                              % F_t
        xx=cat(2,ML_center(F(1:T-1,:))',zeros(r,1));                            % F_{t-1}
        EF1=zeros(r,r); EF=zeros(r,r);  EF1F1=EF;                               % initialize
        for tt=cc+2:T  % ------------------------------------------------------ % \sum_{t=1}^T
            EF1=EF1+yy(:,tt)*xx(:,tt-1)'+PtTm(:,:,tt);                          % E(F_t F_{t-1})
            EF1F1=EF1F1+xx(:,tt-1)*xx(:,tt-1)'+PtT(:,:,tt-1);                   % E(F_{t-1} F_{t-1})
            EF=EF+yy(:,tt)*yy(:,tt)'+PtT(:,:,tt);                               % E(F_t F_t)
        end
        Q = (EF + EF1F1 - 2*EF1') / (T-cc);
    end
    
    %%% ============================== %%%
    %%% 3: Covariance Prediction Error %%%
    %%% ============================== %%%
    yy=y';
    xx=xitT';
    if sr==1
        R=zeros(N,N);                                                           % preallocates
        for tt=cc+1:T  % ------------------------------------------------------ % \sum_{t=1}^T
            eta=yy(:,tt)-Lambda*xx(:,tt);                                       % prediction error
            R = R+ eta*eta'+ Lambda*PtT(:,:,tt)*Lambda';                        % E(eta_t eta_t')
        end
        R = diag(diag(R/(T-cc)));
    end
%     disp([Q R])
    
    %%%% ======================== %%%%
    %%%% ====                ==== %%%%
    %%%% ====    E - Step    ==== %%%%
    %%%% ====                ==== %%%%
    %%%% ======================== %%%%            
    [xitt,xittm,Ptt,Pttm,loglik]=ML_KalmanFilter2(F00,P00,y,A,Lambda,R,Q);  % Kalman Filter
    [xitT,PtT,PtTm]=ML_KalmanSmoother2(A,xitt',xittm',Ptt,Pttm,Lambda,R);   % Kalman SMoother           
    F=xitT;                                                          % jj-th step factors            
    F00=xitT(1,:)'; P00=PtT(:,:,1);                                         % initial conditions    

    %%%% ================================= %%%%
    %%%% ====                         ==== %%%%
    %%%% ====    Check convergence    ==== %%%%
    %%%% ====                         ==== %%%%
    %%%% ================================= %%%%
    delta_loglik = abs(loglik - loglik1);                                   % |logL(t) - logL(t-1)|
    avg_loglik = (abs(loglik) + abs(loglik1) + 10^(-3))/2;                  % avg = (|logL(t)| + |logL(t-1)|)/2
    LogLik(jj+1)=loglik;
    disp([num2str(jj) ': ' num2str(loglik) '  ' ...
        num2str(delta_loglik / avg_loglik)])
    if (delta_loglik / avg_loglik) < tresh; break; end                      % convergence if |f(t) - f(t-1)| / avg < threshold
    if jj>1; if loglik - loglik1<0; disp('decreased'); break; end; end;
    loglik1=loglik;      
end


LLLM.Lambda=Lambda;
LLLM.xitT=xitT;
LLLM.PtT=PtT;
LLLM.PtTm=PtTm;
LLLM.xitt=xitt;
LLLM.xittm=xittm;
LLLM.Ptt=Ptt;
LLLM.Pttm=Pttm;
LLLM.R=R;
LLLM.Q=Q;
