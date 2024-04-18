% BLL_VECM - Estimate a VECM with Reduced Rank Regression as in Johansen
% 
% [AA, alpha, beta, Phi, Theta, u, Gamma, A]=BLL_VECM(x,p,r,s,K)
%  
%   VECM: (1-L)xt = K + alpha*beta'x{t-1} + Gamma(L)(1-L)x{t-1} + ut
%    VAR: A(L)xt=ut
%     MA: xt=Phi(L) ut
% 
%  Inputs:
%  x-
%  p - number of lags in A(L)
%  r - number of cointegration relations
%  s - number of lags to compute for Phi(L)
%  K - constant (yes/no)
% 
%  Outputs:
%     AA - Parameters of A(L) in a Single matrix
%  alpha - 
%   beta - cointegrating vector
%    Phi - Phi(L), Moving average representation
%  Theta - Long Run Theoretical IRF, aka Phi(1)
%      u - ut,residuals
%  Gamma - Gamma(L), Autoregressive parameters first difference representation
%      A - A(L)
%
% Written by Matteo Luciani (matteo.luciani@frb.gov)
% 
% Reference: Barigozzi, M, Lippi, M. and Luciani M. (2016) "Non-Stationary 
%            Dynamic Factor Models for Large Datasets", FEDS 2016-024
%  

function [AA, alpha, beta, Phi, Theta, u, Gamma, A]=BLL_VECM(x,p,r,s,K)

    
[T, n]=size(x);    
if nargin<5; K=0; end                                                       % no input no constant in the VECM

%%% -------------------------------- %%%
%%% Estimate Beta as a RR regression %%%
%%% -------------------------------- %%%
xx=ML_diff(x);                                                              % First Difference
Z0=xx(p:T-1,:);                                                             % \Delta x_t
Z1=x(1:T-p,:);                                                              % x_{t-1}
Z2=ML_lag(xx,p-1);                                                          % \Delta x_{t-i}, i=1:p-1
t=length(Z0);

if p==1;        
    R0=Z0; R1=Z1;                                                           % Residuals
else
    M02=Z0'*Z2/t; M12=Z1'*Z2/t; M22=Z2'*Z2/t;                               % Build Moment Matrices    
    R0=Z0-Z2*inv(M22)'*M02'; R1=Z1-Z2*inv(M22)'*M12';                       % Residuals
end
S00=R0'*R0/t; S01=R0'*R1/t; S10=R1'*R0/t; S11=R1'*R1/t;                     % S Matrices
SS=S11^(-.5)*S10*S00^(-1)*S01*S11^(-.5);
[U, l]=svd(SS); [l, jj]=sort(diag(l)); l=flipud(l); jj=flipud(jj);          % eigenvalue eigenvectors
V=S11^(-.5)*U;                                                              % normalization 
beta=V(:,jj(1:r));           

%%% ----------------------------------------- %%%
%%% Estimate the rest of the parametrs by OLS %%%
%%% ----------------------------------------- %%%
PSI=NaN(n*(p-1)+r+K,n);
Gamma=NaN(n,n,p-1);
u=NaN(T-p,n);
z=x*beta;                                                                   % Error Correction term

for ii=1:n; [PSI(:,ii), u(:,ii)]=ML_ols(xx(p:end,ii),[z(p:end-1,:) Z2],K); end

%%% reshpae the estimated coefficients %%%
if K>1; kk=K-1; else kk=K;end
for pp=1:p-1; Gamma(:,:,pp)=PSI(kk+r+pp:p-1:end,:)'; end 
alpha=PSI(kk+1:kk+r,:)';
PI=alpha*beta';

%%% --------------------- %%%
%%% Implied VAR in Levels %%%
%%% --------------------- %%%
if p==1;
    A(:,:,1) = PI + eye(n);
else
    A(:,:,1) = PI + eye(n) + Gamma(:,:,1);
    if p==2;
        A(:,:,2) = -Gamma(:,:,1);
    else
        for pp=2:p-1;
            A(:,:,pp) = Gamma(:,:,pp) - Gamma(:,:,pp-1);
        end
        A(:,:,p) = -Gamma(:,:,p-1);
    end
end

AA=NaN(n*p,n); for pp=1:p; AA(1+pp-1:p:end,:)=A(:,:,pp)'; end               % VAR coefficient store in a different Way

if kk>0; AA=[PSI(1:kk,:); AA]; end

%%% ----------------- %%%
%%% Impulse Responses %%%
%%% ----------------- %%%
if nargin <5;s=20; end
Phi(:,:,1) = eye(n);
for jj=2:s;
    Phi(:,:,jj) = 0;
    for ii = 1:min(jj-1,p);        
        temp3=A(:,:,ii)*Phi(:,:,jj-ii);        
        Phi(:,:,jj)=Phi(:,:,jj)+temp3;
    end;
end;

%%% --------------------------------------- %%%
%%% Theoretical Long-Run IRF -- C(1) matrix %%%
%%% --------------------------------------- %%%
alphaC=ML_OrthComp(alpha);
betaC=ML_OrthComp(beta);
if p==1;
    Theta=betaC*inv(alphaC'*betaC)*alphaC';
    Gamma=[];
else
    Theta=betaC*inv(alphaC'*(eye(n)-sum(Gamma,3))*betaC)*alphaC';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gammaC=ML_OrthComp(gamma)

[n, r] = size(gamma);
[P,L] = eig(gamma*gamma');
gammaC = P(:,1:n-r);
