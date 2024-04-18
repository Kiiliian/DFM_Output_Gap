function ZZ=BL_Band(x,xhat,confidence)

if nargin <3; confidence = 1; end

if confidence ==5; alpha=1.96;
elseif confidence==1; alpha=2.576;
elseif confidence==10; alpha=1.645;
elseif confidence==16; alpha=1.405;
elseif confidence==32; alpha=0.9945;
end
    
Sigma=std(xhat,[],2);
ZZ=[x-alpha*Sigma x+alpha*Sigma];