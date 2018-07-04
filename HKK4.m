%in the function u is integral variable
%write the equation of our own
%rho=1, because we do dimensionless through rho
% Scattering problem: Seismic dispersion and attenuation in saturated porous rock with aligned slit cracks
% this is a function usd in scatteringproblem, it will be used combining
% with HKK
% This code is a low frequency limit ofthe H in Song (2017)in equation (60)
% All the parameters in thsi function is same as the parameters in
% scatteringproblem
% v.0.1, 18/05/2018, Boye Fu & Boris Gurevich, Curtin University and CRGC
function HKK4 = HKK4(u)
global k0 k1 k2 k3 visc perm b bdim I Kg mug Kf porosity alpha Kdry;
global mudry taudry L Kstar M HBiot Mdim HBiotdim Ldim rhodrydim rhofdim oo;
c0=sqrt(Ldim/rhodrydim);%the dimensionless p-wave velcocity for dry rock
k0=oo/sqrt(Ldim);%the dimensionless P-wave number
c1=sqrt(HBiotdim);%the dimensionless Fast P-wave velocity
k1=oo/c1;%The dimensionless Fast P-wave number
k2=sqrt(I*oo*bdim*HBiotdim/(Ldim*Mdim));%The dimensionless Slow P-wave number
k3=oo;%the demensionless S-wave nuber(the demensionless is relation to shear modulus, therefore, the demensionless shear modulus is 1)
kk=abs(k2);%the abslute value of slow P-wave number
X2=-HBiotdim/(alpha*Mdim);
X3=I*(rhofdim.*oo./bdim);
cigma1=(HBiotdim-alpha*Mdim)/2;
cigma2=Ldim/(2*alpha);
cigma3=(HBiotdim-alpha*Mdim)/2;
cigma4=Mdim*Ldim/(2*HBiotdim*rhofdim);
E=k3^2*(1*(1-2*cigma2)/cigma4-(1+X2)*(1-2*cigma1)/cigma3);
q1=sqrt(u.^2-k1.^2);
q2=sqrt(u.^2-k2.^2);
q3=sqrt(u.^2-k3^2);
%%
H1=4/E;
H2=0;
H3=(0-X3).*((-(cigma2+cigma4).*k2.^2)./q2);
H4=0;
HH=H1.*H3;
HKK4=(besselj(1,u).^2).*HH;%*u/z;%the u/z is because of dimensionless
end