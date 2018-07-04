%%
% Scattering problem: Seismic dispersion and attenuation in saturated porous rock with aligned slit cracks
% this is a function usd in scatteringproblem, it can calculate the f(w) in
% Song (2017) in equation (A13)
% v.0.1, 18/05/2018, Boye Fu & Boris Gurevich, Curtin University and CRGC
%%
function g = gff1(z)
global I kk Kg mug Kf porosity alpha Kdry mudry taudry L Kstar M HBiot;
global Mdim HBiotdim Ldim;
T=besselj(1,z)/z;
p=-1*(HBiotdim-alpha*Mdim);
g=p*T;
end