% BL_StateSpaceLLLM

% Yt = C St + et    et ~ N(0,R)
% St = A St1 + ut   ut ~ N(0,Q)

function SS=BL_StateSpaceLLLM(Y,Q,nu)

if nargin==2; nu=0.95; end
n=size(Y,2);
tau=(1-nu)^2;
SS.A=1;
SS.C=ones(n,1) ;
SS.Q=Q;
SS.R=cov(Y)/length(Y);                                                           

SS.p0=SS.Q/tau; 
SS.f0=mean(ML_diff(mean(Y,2))); 
