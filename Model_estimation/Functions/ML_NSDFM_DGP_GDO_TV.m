% ML_NSDFM_DGP_GDO_TV - Generate data according to I(1) DFM

function [Z1s,xitTs]=ML_NSDFM_DGP_GDO_TV(T,N,q,s,p,EM,eta,type,TREND,block)

A=EM.A; mu=EM.mu; Lambda=EM.Lambda; Q=EM.Q; xitT=EM.xitT; PtT=EM.PtT;       % Maximum Likelihood Estimates
s1=s+1;
r=q*(s1);                                                                   % number of "static" factors
pq=q*p;
pqs=max(pq,r);
ps=max(p,s+1);

isxi=find(type==1)+pqs;                                                     % identifies states that are idiosyncratic components
istrend=find(mu~=0); istrend(1:q)=[];                                       % identifies states with linear trend
isc2=isxi(~ismember(isxi,istrend));                                         % identifies states that are idio RWnLT
isc3=istrend(~ismember(istrend,isxi));                                     	% identifies states that are linear trends
isc4=isxi(ismember(isxi,istrend));                                          % identifies states that are idio RWwLT
J2=find(sum(Lambda(:,isc2),2)==1); nJ2=ismember(1:N,J2);                    % identifies variables that have idio RWnLT 
J3=find(sum(Lambda(:,isc3),2)==1); nJ3=ismember(1:N,J3);                    % identifies variables that have idio WNwLT 
J4=find(sum(Lambda(:,isc4),2)==1); nJ4=ismember(1:N,J4);                    % identifies variables that have idio RWwLT 
isTV=find(type==10)+pqs;                                                    % TV deterministic components
TT=(1:T)';                                                                  % time trend
xitTs=NaN*xitT;                                                             % simulated states

%%% ============================================== %%%
%%% === -------------------------------------- === %%%
%%% === -- Simulate Deterministic Component -- === %%%
%%% === -------------------------------------- === %%%
%%% ============================================== %%%
LT=zeros(T,N);                                                              % preallocates
for ii=1:length(isTV); jj=isTV(ii); % ------------------------------------- % TV slopes or means
    id{ii}=find(Lambda(:,jj));
    if TREND(ii)==1;
        zeta=ML_center(randn(T,1)*sqrt(Q(jj+1,jj+1)));                      % simulated shocks
        Ds(1,1)=xitT(1,jj+1)+randn*sqrt(PtT(jj+1,jj+1,1));                  % initial condition
        for tt=2:T; Ds(tt,1)=Ds(tt-1,1)+zeta(tt,:); end                     % generate TV slope
        LT(:,id{ii})=cumsum(Ds)*ones(1,length(id{ii}));                     % TV Trend
        xitTs(:,jj)=cumsum(Ds); xitTs(:,jj+1)=Ds;                           % simulated states       
    else
        zeta=ML_center(randn(T,1)*sqrt(Q(jj,jj)));                          % simulated shocks
        Ds(1,1)=xitT(1,jj)+randn*sqrt(PtT(jj,jj,1));                        % initial condition
        for tt=2:T; Ds(tt,1)=Ds(tt-1,1)+zeta(tt,:); end                     % generate TV mean
        LT(:,id{ii})=Ds*ones(1,length(id{ii}));                             % TV Mean
        xitTs(:,jj)=Ds;                                                     % simulated states
    end                           
end
LT(:,nJ3)=xitT(:,isc3);                                                     % linear trend for WNwLT
LT(:,nJ4)=TT*mu(isc4)';                                                     % linear trend for RWwLT
xitTs(:,isc3)=xitT(:,isc3);                                                 % simulated states

%%% ======================================= %%%
%%% === ------------------------------- === %%%
%%% === -- Simulate Common Component -- === %%%
%%% === ------------------------------- === %%%
%%% ======================================= %%%
Fs=zeros(T+1,q);                                                            % preallocates
us=randn(T+ps,q)*chol(Q(1:q,1:q));                                          % simulated common shocks

for jj=1:p; Fs(jj,:)=xitT(jj,1:q)+randn(1,q)*chol(PtT(1:q,1:q,jj)); end     % initial condition common factors

for tt=p+1:T+ps;                                                            % Generate Common Factors
    Ftemp=[]; for pp=1:p; Ftemp=cat(1,Ftemp,Fs(tt-pp,:)'); end              % -----------------------
    Fs(tt,:)=mu(1:q)'+(A(1:q,1:pq)*Ftemp)'+us(tt,:);                        % -----------------------
end                                                                         % -----------------------
Fs1=[]; for ss=0:s; Fs1=cat(2,Fs1,Fs(ps-ss+1:end-ss,:));end                 % Lagged factors (s+1) lags for chi
chis=Fs1*Lambda(:,1:r)';                                                    % common component

FsT=[]; for ss=0:ps-1; FsT=cat(2,FsT,Fs(ps-ss+1:end-ss,:));end              % Lagged factors max(p,s+1) lags 
xitTs(:,1:pqs)=FsT;                                                         % simulated states

%%% ============================================== %%%
%%% === -------------------------------------- === %%%
%%% === -- Simulate Idiosyncratic Component -- === %%%
%%% === -------------------------------------- === %%%
%%% ============================================== %%%
EE=zeros(T,N);  Es=EE;                                                      % preallocates
EE(:,nJ2)=xitT(:,isc2);                                                     % idio for RWnLT
EE(:,nJ4)=xitT(:,isc4)-LT(:,nJ4);                                           % idio for RWwLT

J1=[isc2 ; isc4];
J2=[find(nJ2') ; find(nJ4')];
                                            
for jj=1:length(J1); ii=J1(jj); kk=J2(jj);                                  % non stationary part
    v=sqrt(Q(ii,ii))*randn(T+1,1);                                          % idiosyncratic shock   
    Es(1,kk)=EE(1,kk)+sqrt(PtT(ii,ii,1))*randn;                             % first observation
    for tt=2:T; Es(tt,kk)=Es(tt-1,kk)+v(tt); end                            % Recursion for RW
end                                                                         % -------------------

xitTs(:,isc2)=EE(:,nJ2); xitTs(:,isc4)=EE(:,nJ4)+LT(:,nJ4);                 % simulated states


%%% ======================================= %%%
%%% === ------------------------------- === %%%
%%% === -- Simulate stationary part  -- === %%%
%%% === ------------------------------- === %%%
%%% ======================================= %%%
I0=2;
if I0==1; Sigma=cov(eta); etastar=randn(T,N)*chol(Sigma);                   % simulate as normal
else etastar=stationary_bootstrap(eta,block);                               % statiomnary bootstrap
end

%%% ============================= %%%
%%% ===                       === %%%
%%% ===    Simulated Data     === %%%
%%% ===                       === %%%
%%% ============================= %%%

Z1s=LT+chis+Es+etastar;                                                      % simulated data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bsdata, indices]=stationary_bootstrap(data,w)
% Implements the stationary bootstrap for bootstrapping stationary, dependent series
% 
% USAGE:
%   [BSDATA, INDICES] = stationary_bootstrap(DATA,B,W)
% 
% INPUTS:
%   DATA   - T by 1 vector of data to be bootstrapped
%   B      - Number of bootstraps
%   W      - Average block length. P, the probability of starting a new block is defined P=1/W
% 
% OUTPUTS:
%   BSDATA  - T by B matrix of bootstrapped data
%   INDICES - T by B matrix of locations of the original BSDATA=DATA(indexes);
% 
% COMMENTS:
%   To generate bootstrap sequences for other uses, such as bootstrapping vector processes,
%   set DATA to (1:N)'.  
%
% See also block_bootstrap

% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 2    Date: 12/31/2001



t=size(data,1);
p=1/w;                                  % Define the probability of a new block
        % Large p ==> large bias small variance
indices=zeros(t,1);                     % Set up the bsdata and indices
indices(1,:)=ceil(t*rand(1,1));         % Initial positions

% Set up the random numbers
select=rand(t,1)<p;
indices(select)=ceil(rand(1,sum(sum(select)))*t);
for i=2:t % Determine whether we stay (rand>p) or move to a new starting value (rand<p)
    indices(i,~select(i,:))=indices(i-1,~select(i,:))+1;
end
indices(indices>t) = indices(indices>t)-t;
% The indices make finding the bsdata simple
bsdata=data(indices,:);
