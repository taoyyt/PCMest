function J = ZonePcmCostFunction(x,u,e,data)
p=data.PredictionHorizon;

U=u(1:p,data.MVIndex(1))/850.;

X1=x(2:p+1,1);
X=(X1-24-273.15)/10;
X4=max(X,0);
X3=x(2:p+1,3); %PCM temperature

Ta=data.MDIndex(2);

if (X1>297.15) &(X3>=Ta) % cooling or with ambient air
    J=sum(sum(X4.^2));
    
elseif (X1<297.15)&(X3>=Ta) % regeneration
    J=sum(sum(((X3-13-273.15)/10).^2));
    
else % cooling with PCM or no air supply
    J=10*sum(sum(X4.^2))+sum(sum(U.^2));
    
end
