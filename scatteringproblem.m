%%
% Scattering problem: Seismic dispersion and attenuation in saturated porous rock with aligned slit cracks
% The code is got according to Galvin's thesis Galvin, R. (2007), Elastic wave attenuation, dispersion and anisotropy in fractured porous media. Doctoral dissertation, Curtin University.
% The detail of the parameters of the code cn be got from the paper:Song, Y., Hu, H., & Rudnicki, J. W. (2017b), Normal compression wave scattering by a permeable crack in a fluid-saturated poroelastic solid. Acta Mechanica Sinica, 33(2), 356-367.
% v.0.1, 18/05/2018, Boye Fu & Boris Gurevich, Curtin University and CRGC
%
%   The function contain some self-defining function
%   lambdaK_generalf8: Calculate the KB in equation A13 in Song (2017)
%   gff1:calculate the right hand of A13 in Song (2007)
%   lambdaK_generalf:Calculate the R(x,y)T(y) in equation (29) in Galvin
%   (2009):Galvin, R. J., & Gurevich, B. (2009). Effective properties of a poroelastic medium containing a distribution of aligned cracks. Journal of Geophysical Research: Solid Earth, 114(B7).
%   gff: calculate the -p0S(x) in equation (29)in Galvin (2009)
%   Unit of The elastic moduli in the code is Pa, and we do  normalization by the shear modulus of dry rock
%   the unit of density in the code are kg/m3
%%
clear all;
t0=cputime;
global kk_e h tau L_e Ldim_e mu_e k0 k1 k2 k3 kk k2_g visc perm b bdim I;
global rrr oo Kg mug Kf porosity alpha Kdry mudry taudry L Kstar M HBiot Kfdim Kgdim;
global Mdim HBiotdim Ldim rhofdim rhodrydim A1_g A2_g q1 q2 q3;%the modulus for dimensionless 
%%
%Some basic parameters we used in this paper
I=sqrt(-1);% the  imaginary unit 
n0=0.01;% the slit crack density
Kg=37*10^9;%Bulk modulus of rock matrix
mug=44*10^9;%shear modulus of rock matrix
Kf=2.25*10^9;%bulk modulus of pore fluid

porosity=0.1;%porosity
visc=1e-3;%viscous cofficient
perm=1e-13;%permeability
rhof=1e3;%the density of pore fluid
rhos=2.65e3;%density of rock matrix
rhodry=(1-porosity)*rhos;%density of dry rock
rho=(1-porosity)*rhos+porosity*rhof;%density used in Biot equation
b=visc/perm;%the ratio between viscous cofficient and permeability
alpha=1-(1-porosity)^(3/(1-porosity));%the Biot coefficient
Kdry=Kg*(1-alpha);%the dry rock bulk modulus
mudry=mug*(1-alpha);%the dry rock shear modulus
bdim=b/sqrt(rho*mudry);%%I a middle varible-the b for dimensionless
L=Kdry+4*mudry/3;%The P-wave modulus for dry rocks;
taudry=sqrt(mudry/L);%?? a middle variable dimensionless The "g" in Galvin (2009)
M=Kg/((1-Kdry/Kg)-porosity*(1-Kg/Kf));% the Biot modulus of the rock The M in Gassmann equation
HBiot=L+alpha^2*M;%Fast P-wave modulus
Mdim=M/mudry;%The M value for dimensionless
HBiotdim=HBiot/mudry;% the value for dimensionless-Biot Fast p wave modulus
Ldim=L/mudry;% the dimensionless P-wave modulus
Kfdim=Kf/mudry;
rhodrydim=rhodry/rho;%the dimensionless dry rock density
rhofdim=rhof/rho;%the dimensionless dry rock density of rock
Kgdim=Kg/mudry;%the dimensionless bulk modulus of rock matrix
%%

omegadim=[10^-6 10^-5.5 10^-5 10^-4.5 0.0001 10^-3.5 0.001 10^-2.5 10^-2.25 0.01 10^-1.5 0.1 10^-0.5 1 ...
10^0.5 10^0.75 10 10^1.25 10^1.5 100 10^2.25 10^2.5 10^2.75 10^3.25 10^3.5 10000];% frequency
%%
% Now we will calculate the velocity and attenuation in different
% frequency, the method have used by Kawahara and Yamashita (1992) in
% equation (17):Kawahara, J., & Yamashita, T. (1992). Scattering of elastic waves by a fracture zone containing randomly distributed cracks. pure and applied geophysics, 139(1), 121-144.
% according to the method of Kawahara and Yamashita (1992),we first should
% caculate the integral function in equation (A13) in Song (2007)
for q=1:length(omegadim)
    %The basic parameters for calculation 
oo=omegadim(q);% the frequency
c0=sqrt(Ldim/rhodrydim);%the dimensionless p-wave velcocity for dry rock
k0=oo/sqrt(Ldim);%the dimensionless P-wave number
c1=sqrt(HBiotdim);%the dimensionless Fast P-wave velocity
k1=oo/c1;%The dimensionless Fast P-wave number
k2=sqrt(I*oo*bdim*HBiotdim/(Ldim*Mdim));%The dimensionless Slow P-wave number
k3=oo;%the demensionless S-wave nuber(the demensionless is relation to shear modulus, therefore, the demensionless shear modulus is 1)
kk=abs(k2);%the abslute value of slow P-wave number
kk2(q)=kk;
omeg(q)=kk^2;% the the abslute value of slow P-wave number sqrt
h=k1;
x1=0;% the low limit of the integration
if kk <= 1
spac=2*kk;
x2=4;
x2=200*spac;
spac=0.1*kk;
else
spac=0.002*kk;
x2=40;
end
x2=max(1000,kk*50);%the high limit of the integration, becuase the integration region is (0,infinity) in equation (A13),we use the value x2 represent infinity
p=-(HBiotdim-alpha*Mdim);%the f in equation A13
NMAX=200;
omk=eye(NMAX);%the number of discreted net, the detail is in Kawahara and Yamashita (1992) in
% equation (17)
% CALCULATE GAUSS-LEGENDRE WEIGHTS, In this paper, we use GAUSS-LEGENDRE
% integration to deal with the integration in equation (A13)
m=(NMAX+1)/2;%%£¿
xm=(x2+x1)/2;%midpoint 
xl=(x2-x1)/2;%midpoint
NN(q)=xm;
MM(q)=xl;
for i=1:m
Z=cos(pi*(i-1/4)/(NMAX+1/2));%1>Z>0
Z1=Z-1;%0>Z1>-1
ZZ(floor(i))=Z;
ZZZ(i)=Z-Z1;
while (Z-Z1) > 3d-14
p1=1;
p2=0;
for j=1:NMAX
p3=p2;
p2=p1;%=((2*(j-1)-1)*Z*p2-(j-2)*p3)/j
p1=((2*j-1)*Z*p2-(j-1)*p3)/j;%p1=((2*j-1)*Z*((2*(j-1)-1)*Z*p2-(j-2)*p2)/j-(j-1)*p2)/j
end
pp=NMAX*(Z*p1-p2)/(Z^2-1);%z=1 pp=inf z=0, NMAX((j-1)(p2-p1)/j)
Z1=Z;
Z=Z1-p1/pp;
end
x(i)=xm-xl*Z;%x1<x<(x2+x1)/2
x(NMAX+1-i)=xm+xl*Z;%x2>x>(x2+x1)/2
w(i)=2*xl/((1-Z^2)*pp^2);%(x2-x1)/2
w(NMAX+1-i)=w(i);
end
for i=1:NMAX
u(i)=x(i);
end
for j=1:NMAX
z(j)=x(j);
end
for i=1:NMAX
for j=1:NMAX
omkk_g(i,j)=omk(i,j)+lambdaK_generalf8(u(i),z(j))*w(i);%the coefficient matrix of equation (A13), a discretization of the integration of KB in equation (A13) in Song (2017)
end
end
for j=1:NMAX
gmat(j)=gff1(z(j));%The discretization of the right hand in equation (A13) in Song (2017)
end
f_g=linsolve(omkk_g.',gmat.');
B0=f_g(1);% The B value in equation (A13) in Song (2017)
%%
%the s1 s2 s3 are the slowness of fast slow O wave and S wave,which is the low frequency limit of the s1 s2 s3 in Song (2017) 
s1(q)=k1/oo;
s2(q)=k2/oo;
s3(q)=k3/oo;
%%
%X1 X2 X3 is the low frequency limit of the X1 X2 X3 in Song (2017) 
X1(q)=0;
X2(q)=-HBiotdim/(alpha*Mdim);
X3(q)=I*(rhofdim*oo/bdim);
%%
%cigma1 cigma2 cigma3 cigma4 is the low frequency limit of the cigma1
%cigma2 cigma3 cigma4 in Song (2017)  in equation (51) and (53)
cigma1(q)=(HBiotdim-alpha*Mdim)/2;
cigma2(q)=Ldim/(2*alpha);
cigma3(q)=(HBiotdim-alpha*Mdim)/2;
cigma4(q)=Mdim*Ldim/(2*HBiotdim*rhofdim);
E(q)=(2*HBiotdim*alpha*Mdim*(alpha-Ldim)/(Mdim*Ldim*alpha*HBiotdim)+2/(alpha*Mdim)*(1-HBiotdim+alpha*Mdim))*k3^2;%The E value in equation (56) in Song (2017), which is related to calculation and velocity
A0=(3.14*(-X2(q))*(0-cigma3(q)*k1^2)*B0/E(q)/(sqrt(-1)));%The amplitude of P wave according to Galvin (2009)  equation (A6) (34) and (35) the value of A0 is related to velocity and attenuation
Breco(q)=B0;
ff(q)=(I^3)/pi*A0;%;f according to Galvin (2009)  equation (A6) (34) and (35) the value of A0 is related to velocity and attenuation
Q1(q)=4*pi*n0*imag(ff(q));%The attenuation in the whole frequency region
v(q)=1/(1+2*pi*n0*real(ff(q)));%The velocity in the whole frequency region
%%
%Now we will calculate the low frequency limit of the velocity and
%attenuation
Nlow=(2-2*HBiotdim+2*alpha*Mdim)/(alpha*Mdim);% a middle variation to calculate the low frequency limit velocity
ffBBB=@(kk1)HKK(kk1);%The intergation to calculate the low frequency limit attenuation
ffBBBB=@(kk2)HKK4(kk2);%The intergation to calculate the low frequency limit attenuation
FBBc(q)=((integral(ffBBB,0,100000)))+imag((integral(ffBBBB,0,100000)))*1i;
QQ=abs(k2^2);% a middle variation to calculate the low frequency limit velocity
FBB=0.27-1/8*log(QQ);% a middle variation to calculate the low frequency limit velocity
B2=-p*1/Nlow*((1-HBiotdim+alpha*Mdim)/(alpha*Mdim)-Nlow+Ldim/alpha*k2^2*FBB);% The low freuency limit of B in equation (56) in Song (2017)
A2=(pi*(-X2(q))*(0-cigma3(q)*k1^2)*B2/E(q)/(sqrt(-1)));% The low freuency limit of A in equation (56) in Song (2017)
fflow3(q)=(I^3)/pi*A2;%the low frequecny limit f in equation (34)in Galvin (2009)
B2c=-p*(FBBc(q));% The low freuency limit of B in equation (56) in Song (2017)
A2c=(3.14*(-X2(q))*(0-cigma3(q)*k1^2)*B2c/E(q)/(sqrt(-1)));% The low freuency limit of A in equation (56) in Song (2017)
fflow3c(q)=(I^3)/pi*A2c;%the f in equation (34)in Galvin (2009)
vlow3(q)=1/(1+2*3.14*n0*real(fflow3(q)));%The low frequency limit of the velocity
Q1loww3c(q)=4*pi*n0*imag(fflow3c(q));%The low frequency limit of the attenuation
%%
%Now we will caculate the high frequency limit of velocity and attenuation
ffhigh2(q)=sqrt(I)*((HBiotdim-alpha*Mdim)^2)/(2*Mdim*Ldim*abs(k2));%the high frequecny limit f in equation (34)in Galvin (2009)
Q1high2(q)=4*pi*n0*imag(ffhigh2(q));%The high frequency limit of the attenuation
vhigh2(q)=1/(1+2*pi*n0*real(ffhigh2(q)));%The high frequency limit of the velocity

%%
%Branching function
% The code is got according to Guo's paper:Guo, J., Germ¨¢n Rubino, J., Barbosa, N. D., Glubokovskikh, S., & Gurevich, B. (2018), Seismic dispersion and attenuation in saturated porous rocks with aligned fractures of finite thickness: Theory and numerical simulations¡ªPart 1: P-wave perpendicular to the fracture plane. Geophysics, 83(1), WA49-WA62.
T=2*(HBiotdim-alpha*Mdim)^2*(2-4*alpha*taudry^2+3*alpha^2*taudry^4)*n0*bdim/(15*taudry^2*(1-taudry^2)^2*HBiotdim*Ldim);
G=2*pi*n0*(HBiotdim-alpha*Mdim)^2*sqrt((1/(HBiotdim*Mdim*Ldim*bdim)));
C=HBiotdim;
C0=HBiotdim/((1+2*3.14*n0*real(fflow3(q)))^2);
tao=((C-C0)/(C*G))^2;
kesisi=(C-C0)^3/(2*C0*C^2*T*G^2);
Csatt(q)=C./(1+(((C-C0)./C0)./(1-kesisi+kesisi.*sqrt(1-1i.*omegadim(q).*tao./(kesisi^2)))));
vsatt(q)=sqrt(real(Csatt(q)./C));% The velocity of Branch function
QBranch(q)=imag(Csatt(q))./real(Csatt(q));% The attenuation of Branch function
%%
% we will combine Galvin's result
% the part is a copy of Galvin's penny-shaped crack
oo=omegadim(q);
c0=sqrt(Ldim/rhodrydim);
k0=oo/sqrt(Ldim);
c1=sqrt(HBiotdim);
k1=oo/c1;
k2=sqrt(I*oo*bdim*HBiotdim/(Ldim*Mdim));
k3=oo;
kk=abs(k2);
omeg(q)=kk^2;
h=k1;
x1=0;
if kk <= 1
spac=2*kk;
x2=4;
x2=200*spac;
spac=0.1*kk;
else
spac=0.002*kk;
x2=40;
end
x2=max(4,kk*10);
p=-1*(HBiotdim-alpha*Mdim);
NMAX=200;
omk=eye(NMAX);
% CALCULATE GAUSS-LEGENDRE WEIGHTS
m=(NMAX+1)/2;
xm=(x2+x1)/2;
xl=(x2-x1)/2;
for i=1:m
Z=cos(pi*(i-1/4)/(NMAX+1/2));
Z1=Z-1;
while (Z-Z1) > 3d-14
p1=1;
p2=0;
for j=1:NMAX
p3=p2;
p2=p1;
p1=((2*j-1)*Z*p2-(j-1)*p3)/j;
end
pp=NMAX*(Z*p1-p2)/(Z^2-1);
Z1=Z;
Z=Z1-p1/pp;
end
x(i)=xm-xl*Z;
x(NMAX+1-i)=xm+xl*Z;
w(i)=2*xl/((1-Z^2)*pp^2);
w(NMAX+1-i)=w(i);
end
for i=1:NMAX
u(i)=x(i);
end
for j=1:NMAX
z(j)=x(j);
end
for i=1:NMAX
for j=1:NMAX
omkk_g(i,j)=omk(i,j)-lambdaK_generalf(u(i),z(j))*w(i);
end
end
for j=1:NMAX
gmat(j)=gff(z(j));
end
f_g=linsolve(omkk_g.',gmat.');
B0=2*f_g(1)/(pi);
A0=-1*(1-alpha*Mdim/HBiotdim)*B0/(2*(1-taudry^2));
kappa=1+2*pi*n0*A0;
vG(q)=1/real(kappa);
QG(q)=-2*imag(kappa)/real(kappa);
end
%%
figure(1); 
plot(log10(kk2.^2),log10(QG),'*');
hold on;
plot(log10(kk2.^2),log10(Q1),'.')
hold on;
plot(log10(kk2.^2),log10(Q1high2),'r')
hold on
plot(log10(kk2.^2),log10(Q1loww3c),'b')
hold on
plot(log10(kk2.^2),log10(QBranch),'k')
legend('Penny-shaped crack','Slit crack','High frequency asymptotic solution','Low frequency asymptotic solution','Branching function')
xlabel ('Frequency')
ylabel ('Log Attenuation')
grid on;
grid on;
figure(2); 
plot(log10(kk2.^2),vG,'*');
hold on;
plot(log10(kk2.^2),v,'.')
hold on;
plot(log10(kk2.^2),vhigh2,'r')
hold on;
plot(log10(kk2.^2),vlow3,'g')
hold on;
plot(log10(kk2.^2),vsatt,'b')
legend('Penny-shaped crack','Slit crack','High frequency asymptotic solution','Low frequency asymptotic solution','Branching function')
xlabel ('Frequency')
ylabel ('Velocity')
grid on;
grid on;
%calculate new value
LQG=log10(-QG);
LQ1=log10(Q1);
LKK2=log10(omegadim);
NN=length(LKK2);
dLQG=zeros(1,NN-1);
dLQ1=zeros(1,NN-1);
for i=1:NN-1;
    dLQG(i)=(LQG(i+1)-LQG(i))./((LKK2(i+1)-LKK2(i)));
    dLQ1(i)=(LQ1(i+1)-LQ1(i))./((LKK2(i+1)-LKK2(i)));
end
figure(17)
plot(LKK2(1:NN-1),dLQG,'*')
hold on
plot(LKK2(1:NN-1),dLQ1,'*')
legend('Penny-shaped crack','Slit crack')
xlabel ('Frequency')
ylabel ('Log Attenuation')
t1=cputime-t0
t1=cputime-t0