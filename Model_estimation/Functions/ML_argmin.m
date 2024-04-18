% ml_argmin - select the position of the minimum value in x

% Written Matteo Luciani (matteoluciani@yahoo.it)

function k=ML_argmin(x,dim)
if nargin == 1; dim =1; end
if dim==2; x=x'; end
n=size(x,2);
k=zeros(1,n);
for jj=1:n;
    [~,k(1,jj)]= min(x(:,jj));
end;


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 %%%%
%%%%   Old Version   %%%%
%%%%                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% [kmax,n]=size(x);
% k=zeros(1,n);
% for jj=1:n;
% %     jj
% %     k(1,jj)=find(x(:,jj)==min(x(:,jj)));
%     for ii=1:kmax;
%         if x(ii,jj)==min(x(:,jj)); k(1,jj)=ii;  end;        
%     end 
% end;