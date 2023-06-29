simtime = 10;
res = sim('DC_motor.slx', simtime);

Td = 1e-4;
x1 = res.x1.Data;
x2 = res.x2.Data;
x3 = res.x3.Data;
y = res.x2.Data;

size_t = size(x2);

x1_ = zeros(size_t);
x2_ = zeros(size_t);
x3_ = zeros(size_t);
Y = zeros(size_t);

for i = 2 : size_t
    Y(i)=y(i)-y(i-1);
    x1_(i)=(x1(i)+x1(i-1))*Td/2;
    x2_(i)=(x2(i)+x2(i-1))*Td/2;
    x3_(i)=(x3(i)+x3(i-1))*Td/2;
end

X = [x1_, x2_, x3_];
K = (X'*X)^(-1)*X'*Y

La_r=1/(K(1))
Ra_r=-K(2)*La_r
Psi_r=-K(3)*La_r