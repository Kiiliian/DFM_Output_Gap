% ML_Spectrum - Computes the spectrum of a time series
%
% [s omega]=ML_Spectrum(x,M,method)
%

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)
% Modified version of spectral.m by Matteo Barigozzi

% spec - bartlett estimation of spectra; 
% m is the lag window size
% the macro is based on crosspec.m

function [s, omega]=ML_Spectrum(x,M,graph,cs)


[T, N]=size(x);
s=NaN(65,N);

if nargin==1; M = floor(sqrt(T)); elseif isempty(M); M = floor(sqrt(T));  end
if nargin<3; graph=0; end
if nargin<4; cs=0; end

if cs==1;
    for nn=1:N
        for jj=1:N;
            [c, omega]=crosspec(x(:,nn),x(:,jj),M);
            s(:,nn,jj)=real(c(1:65));
        end
    end
else
    for nn=1:N
        [c, omega]=crosspec(x(:,nn),x(:,nn),M);
        s(:,nn)=real(c(1:65));
    end
end
omega=omega(1:65);

if graph==1;    
    figure; ML_FigureSize(2);
    if N==1;
        plot(omega,s,'linewidth',1.5); axis tight
    else
        plot(omega,s./repmat(max(s),65,1),'linewidth',1.5); axis tight
    end
    set(gca,'Xtick', omega(1:5:end));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crosspec - estimation of cross-spectra.
%
%            c=crosspec(x,y,m,method)
%
% m is the lag window size;  method is the window type which
% can be 'bartlett', 'hamming', or 'hanning' (other windows
% can be constructed by the user). If method is unspecified
% a bartlett window is used.
% frequencies are calculated at 128 equally
% spaced points between 0 and 2*pi-2*pi/128.

function [c, omega]=crosspec(x,y,m,method)
omega=0:2*pi/128:2*pi-2*pi/128;
if nargin==3,
    method='bartlett';
end
C=crosscov(x,y,m).*feval(method,2*m+1);
c=(C'*exp(-i*(-m:m)'*omega))';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crosscov - 
% usage: y=crosscov(x,z,n)
% is similar to xcov but truncates cross-covariances
% at lag n and lead n. The 'unbiased' version of xcov is used

function y=crosscov(x,z,n)
T=length(x);
if nargin==2,
    n=z;
    b=xcov(x,'unbiased')*T/(T-1);
else
    b=xcov(x,z,'unbiased')*T/(T-1);
end;
c=length(b(:,1));
t=(c+1)/2;
tt=(t+n):-1:(t-n);
y=b(tt,:);
