%E.Drive Lab1 var = 13
%Task 1

m = 2200;
v = 0.013;
z1 = 13;
z2 = 28;
z3 = 17;
z4 = 40;
z5 = 17;
z6 = 33;
J2 = 0.12;
J3 = 0.2773;
J4 = 0.0693;
J5 = 0.3216;
J6 = 0.0907;
J7 = 0.4595;
J8 = 0.1040;
J9 = 0.5570;

d1 = 150/1000;
d2 = 200/1000;
j23prev = d2/d1;
z = 1;
s = 0.01;
nr = 0.95;
nz = 0.9;
nv = 0.6;
nc = 0.97;
jc = 1;


cm = 15*10^6;
cr = 3*10^6;
cc = 12*10^6;

%%
% v = w * s * z/2pi
w5 = v*2*pi/(s*z);
j23 = d2/d1; 
j45 = z2/z1; 
j67 = z4/z3; 
j89 = z6/z5;
wd = d2/d1*z2/z1*z4/z3*z6/z5*w5
nd = 30*wd/pi
Md = m/2*10*v*(1+1/nc)/(wd*nr*nz*nz*nz*nv)
Pd=Md*wd
Pd_kvt =Pd/1000 
%%
%Motor  - 5A80MA6
Pnom = 0.75*1000;
nnom = 930;
Mnom = 7.7;
wnom = pi*nnom/30;
nnomh=1000;
wnomh = pi*nnomh/30;

w2 = wd*d1/d2;
n2 = nd*d1/d2;
M2 = Md*nr*d2/d1;
h = Mnom/(wnomh-wnom);
% h = Mnom/( nnomh-nnom);
%%
d1 = 150/1000;
d2 = 200/1000;
j23prev = d2/d1;
%%
j23new1 = wnomh/(2*w2)*(1+sqrt(1-4*M2*w2/(nr*h*wnomh^2)))
% j12new1 = nnomh/(2*n2)*(1+sqrt(1-4*M2*n2/(nr*h*nnomh^2)))


%j12new2 = wnomh/(2*w2)*(1-sqrt(1-4*M2*w2/(nr*h*wnomh^2)))
d1=d2/j23new1
% wd = d2/d1*z2/z1*z4/z3*z6/z5*w5
% nd = 30*wd/pi
% Md = m/2*10*v*(1+1/nc)/(wd*nr*nz*nz*nz*nv)
% Pd=Md*wd
% Pd_kvt =Pd/1000 
%%
%Task2
%a)
Jd = 0.0033;
J37_= (J3+J4)/j23 + (J5+J6)/(j23*j45)+J7/(j23*j45*j67)
J810_= J8/(j23*j45*j67) + J9/(j23*j45*j67*j89)+m/2*(v/wd)^2
J10_2_= m/2*(v/wd)^2

c2_ = c2
c3_ = c3/(j23*j45*j67)^2
c4_ = c4/(j23*j45*j67*j89)^2

%b
J17_ = Jd + J2 + (J3+J4)/j23 + (J5+J6)/(j23*j45)+J7/(j23*j45*j67)



%c
J810_2= J8/(j23*j45*j67) + J9/(j23*j45*j67*j89)+m*(v/wd)^2

c34_ = (1/c3_+1/c4_)^(-1)

%%
% Task3

