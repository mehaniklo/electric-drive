function mu_m = mu_m_calc(r_2,r1,s_n,lyambda,In,cosfi,Un,m,z_p,w1,Mn)
    a = r1/r_2;
    A = 1-2*a*s_n*(lyambda-1);
    s_m = s_n*(lyambda+sqrt(lyambda^2-A))/A;
    x_ks = sqrt((r_2/s_m)^2-r1^2);
    b = x_ks/((r1+r_2/s_n)^2+x_ks^2);
    x_m = 1/((In*sqrt(1-cosfi^2))/Un-b);
    I2_ = Un/sqrt((r1+r_2/s_m)^2+x_ks^2);
    mu_m = m*z_p*abs(I2_^2)*r_2/(w1*s_m*Mn)-lyambda; 
end