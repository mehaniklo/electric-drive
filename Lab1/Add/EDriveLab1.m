1%E.Drive Lab1 var = 13
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
w5 = v*2*pi/(s*z)
j23 = d2/d1
j45 = z2/z1 
j67 = z4/z3 
j89 = z6/z5
wd = d2/d1*z2/z1*z4/z3*z6/z5*w5

nd = 30*wd/pi
Md = m/2*9.81*v*(1+1/nc)/(wd*nr*nz*nz*nz*nv)
Pd=Md*wd
Pd_kvt =Pd/1000 

ML1=m/2*9.81*v/(wd*nr*nz*nz*nz*nv);
ML2=m/2*9.81*v/(wd*nr*nz*nz*nz*nv*nc);

%%
%Motor  - 5A80MA6
Pnom = 0.75*1000
nnom = 930
Mnom = 7.7
wnom = pi*nnom/30
nnomh=1000
wnomh = pi*nnomh/30

w2 = wd*d1/d2
n2 = nd*d1/d2
M2 = Md*nr*d2/d1
h = Mnom/(wnomh-wnom)
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
%%
%Task2
%a)
Jd = 0.0033;
J37_= (J3+J4)/j23^2 + (J5+J6)/(j23*j45)^2+J7/(j23*j45*j67)^2
J810_= J8/(j23*j45*j67)^2 + J9/(j23*j45*j67*j89)^2+m/2*(v/wd)^2
J10_2_= m/2*(v/wd)^2

c1_ = cm
c2_ = cr
c3_ = cm/(j23*j45*j67)^2
c4_ = cc/(j23*j45*j67*j89)^2

%b
J17_ = Jd + J2 + (J3+J4)/j23^2 + (J5+J6)/(j23*j45)^2+J7/(j23*j45*j67)^2



%c
J810_2_= J8/(j23*j45*j67)^2 + J9/(j23*j45*j67*j89)^2+m*(v/wd)^2

c34_ = (1/c3_+1/c4_)^(-1)
%%
% резонанс
%a)
n = 5;
J1 = Jd;
J2;
J3 = J37_;
J4 = J810_;
J5 = J10_2_;
c1 = cm;
c2 = cr;
c3 = c3_;
c4 = c4_;

A =[0 0 0 0 0  -1/J1   0     0     0;
    0 0 0 0 0   1/J2 -1/J2   0     0;
    0 0 0 0 0    0    1/J3 -1/J3   0;
    0 0 0 0 0    0     0    1/J4 -1/J4;
    0 0 0 0 0    0     0      0   1/J5;
    c1 -c1   0   0   0 0 0 0 0;
    0   c2 -c2   0   0 0 0 0 0;
    0    0  c3 -c3   0 0 0 0 0;
    0    0   0  c4 -c4 0 0 0 0;]
reson5 = eig(A)
%%
%b)3-х массовая
reson3 = eig([0 0 0  -1/J17_             0;
              0 0 0   1/J810_    -1/J810_;
              0 0 0      0        1/J10_2_; 
              c3_ -c3_ 0 0 0;
              0 c4_ -c4_ 0 0])
%a = (J10_2_*c3_*(J17_+J810_)+J17_*c4_*(J810_+J10_2_))/(J17_*J810_*J10_2_);
%b = c3_*c4_*(J17_+J810_+J10_2_)/(J17_*J810_*J10_2_);

%omega3m1 = sqrt(a/2*(1-sqrt(1-4*b/a^2)))
%omega3m2 = sqrt(a/2*(1+sqrt(1-4*b/a^2)))
%c)2-массовая
%omega2m = sqrt(c34_*(J17_+J810_2_)/(J17_*J810_2_))  
reson2 = eig([ 0     0 -1/J17_; ...
               0     0  1/J810_2_; ...
              c34_ -c34_  0])
%% Task3
% 2 mass
Mcf = 0.1; %Mcf/(n1n*j1n)
Mcf1_ =  Mcf;
Mcf2_ = Mcf/(nr*j23)
Mcf7_ = Mcf/(nr*nz*nz*j23*j45*j67)
Mcf10_1_ = Mcf*w5/(nr*nz*nz*nz*nv*wd)
Mcf10_2_ = Mcf10_1_/nc
%2 mass
Mcf1=Mcf1_+Mcf2_+Mcf7_
Mcf2 = Mcf10_1_+Mcf10_2_
J1 = J17_
J2 = J810_2_
K12=c34_
%%
simtime = 20;
res = sim('E_DriveLab1Mass2Simscape.slx', simtime);
% res = sim('E_DriveLab1Mass2.slx', simtime);
%%
u = res.out.signals.values(:,1);
w1 = res.out.signals.values(:,2);
w2 = res.out.signals.values(:,3);
w1_wdf = res.out.signals.values(:,4);
w2_wdf = res.out.signals.values(:,5);
t = res.out.time;


figure;
%tiledlayout(2,1)

%nexttile
plot(t, w1, 'b', 'LineWidth', 2, 'DisplayName','w1')
grid;
xlabel('t');
ylabel("w(t)");
title('Угловая скорость w(t) без учета сухого трения');
%lgd = legend;
%lgd.NumColumns = 1;
%lgd.Location = 'northeast';

%nexttile
hold on
plot(t, w2, 'r', 'LineWidth', 1, 'DisplayName','w2')
%grid;
%xlabel('t');
%ylabel("w(t)");
%title('Угловая скорости w2(t)');
hold off
lgd = legend;
lgd.NumColumns = 1;
lgd.Location = 'northeast';

figure;
%tiledlayout(2,1)

%nexttile
plot(t, w1_wdf, 'b', 'LineWidth', 2, 'DisplayName','w1')
grid;
xlabel('t');
ylabel("w(t)");
title('Угловая скорость w(t) с учетом сухого трения');
%lgd = legend;
%lgd.NumColumns = 1;
%lgd.Location = 'northeast';

%nexttile
hold on
plot(t, w2_wdf, 'r', 'LineWidth', 1, 'DisplayName','w2')
%grid;
%xlabel('t');
%ylabel("w(t)");
%title('Угловая скорости w2(t)');
hold off
lgd = legend;
lgd.NumColumns = 1;
lgd.Location = 'northeast';

figure;
plot(t, u, 'b', 'LineWidth', 1, 'DisplayName','u')
grid;
xlabel('t');
ylabel("u(t)");
title('Внешнее воздействие u(t)');

lgd = legend;
lgd.NumColumns = 1;
lgd.Location = 'northeast';
%% 3mass
Mcf1=Mcf1_+Mcf2_+Mcf7_
Mcf2 = Mcf10_1_
Mcf3 = Mcf10_2_

J1 = J17_

J2 = J810_
J3 = J10_2_


K12=c3_
K23 = c4_

%%
simtime = 20;
% res = sim('E_DriveLab1Mass3.slx', simtime);
res = sim('E_DriveLab1Mass3Simscape.slx', simtime);
%%
u = res.out.signals.values(:,1);
w1 = res.out.signals.values(:,2);
w2 = res.out.signals.values(:,3);
w3 = res.out.signals.values(:,4);
w1_wdf = res.out.signals.values(:,5);
w2_wdf = res.out.signals.values(:,6);
w3_wdf = res.out.signals.values(:,7);
t = res.out.time;


figure;
%tiledlayout(2,1)

%nexttile
plot(t, w1, 'b', 'LineWidth', 1, 'DisplayName','w1')
grid;
xlabel('t');
ylabel("w(t)");
title('Угловая скорость w(t) без учета сухого трения');


%nexttile
hold on
plot(t, w2, 'g', 'LineWidth', 1, 'DisplayName','w2', "LineStyle", '--')

plot(t, w3, 'r', 'LineWidth', 1, 'DisplayName','w3')

hold off
lgd = legend;
lgd.NumColumns = 1;
lgd.Location = 'northeast';

figure;
%tiledlayout(2,1)

%nexttile
plot(t, w1_wdf, 'b', 'LineWidth', 1, 'DisplayName','w1')
grid;
xlabel('t');
ylabel("w(t)");
title('Угловая скорость w(t) с учетом сухого трения');


%nexttile
hold on
plot(t, w2_wdf, 'g', 'LineWidth', 1, 'DisplayName','w2', 'LineStyle', '--')

plot(t, w3_wdf, 'r', 'LineWidth', 1, 'DisplayName','w3')

hold off
lgd = legend;
lgd.NumColumns = 1;
lgd.Location = 'northeast';

figure;
plot(t, u, 'b', 'LineWidth', 1, 'DisplayName','u')
grid;
xlabel('t');
ylabel("u(t)");
title('Внешнее воздействие u(t)');

lgd = legend;
lgd.NumColumns = 1;
lgd.Location = 'northeast';

%% 5 mass
%J1 = J1;
%J2 = J2;
J3 = J37_
J4 = J810_
J5 = J10_2_

K12 = c1_
K23 = c2_
K34 = c3_
K45 = c4_

Mcf1 = Mcf1_;
Mcf2 = Mcf2_;
Mcf3 = Mcf7_;
Mcf4 = Mcf10_1_;
Mcf5 = Mcf10_2_;
%%

simtime = 20;
% res = sim('E_DriveLab1Mass5.slx', simtime);
res = sim('E_DriveLab1Mass5Simscape.slx', simtime);

u = res.out.signals.values(:,1);
w1 = res.out.signals.values(:,2);
w2 = res.out.signals.values(:,3);
w3 = res.out.signals.values(:,4);
w4 = res.out.signals.values(:,5);


w5 = res.out.signals.values(:,6);
w1_wdf = res.out.signals.values(:,7);
w2_wdf = res.out.signals.values(:,8);
w3_wdf = res.out.signals.values(:,9);
w4_wdf = res.out.signals.values(:,10);
w5_wdf = res.out.signals.values(:,11);
t = res.out.time;

figure;
%tiledlayout(2,1)

%nexttile
plot(t, w1, 'b', 'LineWidth', 1, 'DisplayName','w1')
grid;
xlabel('t');
ylabel("w(t)");
title('Угловая скорость w(t) без учета сухого трения');

%nexttile
hold on
plot(t, w2, 'm', 'LineWidth', 1, 'DisplayName','w2', 'LineStyle','-.')
plot(t, w3, 'g', 'LineWidth', 1, 'DisplayName','w3', "LineStyle", '--')
plot(t, w4, 'c', 'LineWidth', 1, 'DisplayName','w4','LineStyle', '--')
plot(t, w5, 'r', 'LineWidth', 1, 'DisplayName','w5','LineStyle', '-')
% "-"	Solid line	
% Sample of solid line
% 
% "--"	Dashed line	
% Sample of dashed line
% 
% ":"	Dotted line	
% Sample of dotted line
% 
% "-."	Dash-dotted line	

hold off
lgd = legend;
lgd.NumColumns = 1;
lgd.Location = 'northeast';

figure;
%tiledlayout(2,1)

%nexttile
plot(t, w1_wdf, 'b', 'LineWidth', 1, 'DisplayName','w1')
grid;
xlabel('t');
ylabel("w(t)");
title('Угловая скорость w(t) с учетом сухого трения');


%nexttile
hold on
plot(t, w2_wdf, 'm', 'LineWidth', 1, 'DisplayName','w2', 'LineStyle','-.')
plot(t, w3_wdf, 'g', 'LineWidth', 1, 'DisplayName','w3', 'LineStyle', '--')
plot(t, w4_wdf, 'c', 'LineWidth', 1, 'DisplayName','w4','LineStyle', '--')
plot(t, w5_wdf, 'r', 'LineWidth', 1, 'DisplayName','w5','LineStyle', '-')

hold off
lgd = legend;
lgd.NumColumns = 1;
lgd.Location = 'northeast';

figure;
plot(t, u, 'b', 'LineWidth', 1, 'DisplayName','u')
grid;
xlabel('t');
ylabel("u(t)");
title('Внешнее воздействие u(t)');

lgd = legend;
lgd.NumColumns = 1;
lgd.Location = 'northeast';

%% Additional tasks
% Дополнительное задание 1
% Идентификация параметров
% Исходные данные: файл, содержащий записи момента и скорости первой 
% массы с периодом дискретизации 1 мс.
% Задания: 
% а) построить частотную характеристику по экспериментальным данным; 
% б) определить параметры двухмассовой механической системы; 
% в) собрать модель с идентифицированными параметрами; 
% г) промоделировать с входным воздействием из исходного файла; 
% д) посчитать среднеквадратичную ошибку между скоростью первой массы с 
% модели и скоростью первой массы из исходного файла.

u = sim_res.out_var.signals.values(:,1);
w1 = sim_res.out_var.signals.values(:,2);
t = sim_res.out_var.time;
%%
figure;
ax1 = subplot(2,1,1);
plot(t, u)
grid on
xlabel("t, с")
ylabel("М, Н*м")
ax2 = subplot(2,1,2);
plot(t, w1)
grid on
xlabel("t, с")
ylabel("w1(t), рад/с")
linkaxes([ax1,ax2],'x')
%%
[txy, f] = tfestimate(u, w1, [], [], [], 1000)
%%
om_1 = 2*pi*f;
L_1=20*log10(abs(txy));
%%
figure();
plot(om_1,L_1);
xlim([0,3000]);
xlabel('\omega_1, rad/s')
ylabel('L, dB')
ax = gca;
set(ax,'Xscale', 'log')
grid on
hold all
%% Подбор параметров модели
Lob = -36.3;
%Lo%b = -25.0007;
om1 = 278.61214;
%om2 = 718.862;
om2 = 737;
Kw = 10^(Lob/20);
T1 = 1/om1;
T2 = 1/om2;
%Попробуем подобрать коэф. демфирования
ksi1 = 0.00321654;%при повышении амплитуда 1 р.ч. увеличвается вверх
ksi2 = 0.0855;%при повышении амплитуда 2 р.ч уменьшается вниз
% ksi1 =0.05;
% ksi2 = 0.05;

%%
%Передаточная функция
s = tf('s')

W_ob = Kw*(T1^2*s^2 + 2*T1*ksi1*s +1)/s/(T2^2*s^2 + 2*T2*ksi2*s +1);

bodemag(W_ob,'r')
grid on
%%
J1=T2^2/(T1^2*Kw)
J2=1/Kw-J1
K12=J2/T1^2
B12=2*ksi2*T2*K12
%% E_DriveLab1AddTask modeling
simtime = 250;
% res = sim('E_DriveLab1Mass5.slx', simtime);
res = sim('E_DriveLab1AddTask.slx', simtime);
%%
u_counted = res.out.signals.values(:,1);
w1_counted = res.out.signals.values(:,2);
t_counted = res.out.time;
%%
figure;
ax1 = subplot(2,1,1);
plot(t, w1)
title("Графики моделирования исходных даных и вычисленной модели")
grid on
xlabel("t, с")
ylabel("w1(t), рад/с")
ax2 = subplot(2,1,2);
plot(t_counted, w1_counted)
grid on
xlabel("t, с")
ylabel("w1(t) counted, рад/с")
linkaxes([ax1,ax2],'x')
%%
e_mse = mse(w1, w1_counted)
%% add Task2
% Исходные данные: файл, содержащий записи момента и скорости первой 
% массы с периодом дискретизации 1 мс и значение коэффициента жесткости 
% K12.
% Задания: 
% а) определить параметры двухмассовой механической системы; 
% б) собрать модель с идентифицированными параметрами; 
% в) промоделировать с входным воздействием из исходного файла; 
% г) посчитать среднеквадратичную ошибку между скоростью первой массы с 
% модели и скоростью первой массы из исходного файла
u = sim_res.out_var.signals.values(:,1);
w1 = sim_res.out_var.signals.values(:,2);
t = sim_res.out_var.time;
%%
figure;
ax1 = subplot(2,1,1);
plot(t, u)
grid on
ax2 = subplot(2,1,2);
plot(t, w1)
grid on
linkaxes([ax1,ax2],'x')
%%
w1_osc = detrend(w1, 1, u == -8);

figure;
ax1 = subplot(2,1,1);
plot(t, u)
grid on
ax2 = subplot(2,1,2);
plot(t, w1_osc)
grid on
linkaxes([ax1,ax2],'x')
%%
Kw = 0.0656555/1.089/8

T2=(2.155-2.03)/(2*pi)
ksi2=-log(0.00559427/0.00831248)/(2*pi)
%%
B12 = ksi2*T2*2*K12

syms x;
J_roots = solve(Kw*x^2 - x + T2^2*K12 == 0, x);
double(J_roots)
J2 = double(J_roots(2))
T1=sqrt(J2/K12)
ksi1=B12/(2*K12*T1)

J1=1/Kw-J2
%% E_DriveLab1AddTask2 modeling
simtime = 16;
res = sim('E_DriveLab1AddTask.slx', simtime);
%%
u_counted = res.out.signals.values(:,1);
w1_counted = res.out.signals.values(:,2);
t_counted = res.out.time;
%%
figure
plot(t, w1)
%%
e_mse = mse(w1, w_counted)