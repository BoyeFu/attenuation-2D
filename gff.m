%%
%%
% Scattering problem: Seismic dispersion and attenuation in saturated porous rock with aligned slit cracks
% this is a function usd in scatteringproblem, it can calculate the pS in
% GAlvin (2009) in equation (29)
% v.0.1, 18/05/2018, Boye Fu & Boris Gurevich, Curtin University and CRGC
%%
function g = gff(z)
global I kk Kg mug Kf porosity alpha Kdry mudry taudry L Kstar M HBiot;
global Mdim HBiotdim Ldim;
T=(sin(z)-z*cos(z))/z^3;
p=-1*(HBiotdim-alpha*Mdim);
g=p*T;
end