% ML_efactors2 - Static Factors Estimation by Principal Components
%
% [fhat,lambda,chat,ehat]=ML_efactors(y,r,jj);
%   Inputs:
%       y is data matrix
%       r is the estimated number of true factors
%       jj is the normalization used in the estimation:
%       jj=1 use TxT, then F'*F/T=I
%       jj=2 us NxN, then lambda'*lamnda/N=I
%   Outputs:
%       fhat is T by r matrix of estimated factors
%       lambda N by r is the estimated loading matrix
%       chat=F*Lambda' T by N is is the estimated common component
%       ehat is T by N matrix of ideosyncratic components
%       V is r by r diagonal matrix containing the first r eigenvalues
%

% Written by Matteo Luciani (matteoluciani@yahoo.it)
% this is a modified version of the codes available on Serena Ng webpage

function [fhat,lambda,chat,ehat,V]=ML_efactors2(y,r,jj)

[T, N]=size(y);
if nargin==2; if T>N; jj=2; else; jj=1; end; end
% w=inv(T*N);
opt.disp = 0;

if jj==1;
%     G=w*y*y'; G=G/G(1,1);
    G=cov(y');
    [dfhat0,eigval] = eigs(G, r,'LM',opt);
    fhat=dfhat0(:,1:r)*sqrt(T); 
    lambda=(fhat'*y/T)';  
    chat=fhat*lambda'; 
    ehat=y-chat;
    V=eigval(1:r,1:r);
elseif jj==2; 
%     G=w*y'*y; G=G/G(1,1);
    G=cov(y);
    [W,M] = eigs(G, r,'LM',opt);
    lambda=W*sqrt(N); 
    fhat=y*W/sqrt(N); 
    chat=fhat*lambda';
    ehat=y-chat;
    V=M(1:r,1:r);
end;


ii=find(sign(fhat(1,:))<0); fhat(:,ii)=-fhat(:,ii); lambda(:,ii)=-lambda(:,ii); % Normalization





