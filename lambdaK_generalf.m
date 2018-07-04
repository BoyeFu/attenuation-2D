%%
% Scattering problem: Seismic dispersion and attenuation in saturated porous rock with aligned slit cracks
% this is a function usd in scatteringproblem, it can calculate the R(x,y)T(y) in
% Galvin (2009) in equation (29)
% v.0.1, 18/05/2018, Boye Fu & Boris Gurevich, Curtin University and CRGC
%%
function lambdaK_general = lambdaK_generalf(u,z)
global k0 k1 k2 k3 visc perm b bdim I oo Kg mug Kf porosity alpha Kdry;
global mudry taudry L Kstar M HBiot Mdim HBiotdim Ldim rhodrydim;
if z==u
K=(-1/pi)*(1-sin(z+u)/(z+u));% The function k in Galvin's Paper
else
K=(-1/pi)*(sin(z-u)/(z-u)-sin(z+u)/(z+u));
end
k0=oo/sqrt(Ldim);
c1=sqrt(HBiotdim);
k1=oo/c1;
k2=sqrt(I*oo*bdim*HBiotdim/(Ldim*Mdim));
k3=oo;
if u<=k1
q1=-I*sqrt(k1^2-u^2);
else q1=sqrt(u^2-k1^2);
end
q2=-I*sqrt(k2^2-u^2);
if u<=k3
q3=-I*sqrt(k3^2-u^2);
else q3=sqrt(u^2-k3^2);
end
H1=1+alpha*Mdim*oo^2/(HBiotdim*(2*u^2-oo^2));%the factor multiplied with T1 and T2
H2=(2*u^2-oo^2)^2-4*u^2*q3*q1-(k1^2-k0^2)*(2*u^2-oo^2)/(alpha*taudry^2);%the numerater of T2
H3=-2*(1-taudry^2)*oo^2*q1*u;%The denominator of T2
H4=2*u^2*alpha^2*Mdim-HBiotdim*bdim*alpha*I*oo;%the numerater of T1
H5=(2*u^2)^2-4*u^2*q3*q2-1*k2^2*(2*u^2-oo^2)/(alpha*taudry^2);
%H5=(2*u^2-oo^2)^2-4*u^2*q3*q2-1*k2^2*(2*u^2-oo^2)/(alpha*taudry^2);
H6=2*HBiotdim*(Ldim-1)*q2*u*(-1)*k2^2;%The factor of denominater of T1
H7=2*u^2-oo^2*(1-alpha*Mdim/HBiotdim);%The factor of denominater of T1
H=H1*(H2/H3+H4*H5/(H6*H7))-1;
lambdaK_general=K*H*u/z;%the u/z is because of dimensionless
end