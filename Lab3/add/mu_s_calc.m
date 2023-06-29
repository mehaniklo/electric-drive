function mu_s = mu_s_calc(h,m,z_p,Un,r2_,w1,Mn,r1,x_s1,x_s2,ks)
    k_r = h*(sinh(2*h)+sin(2*h))/(cosh(2*h)-cos(2*h));
    k_x = 3/(2*h)*(sinh(2*h)-sin(2*h))/(cosh(2*h)-cos(2*h));
    mu_s = m*z_p*Un^2*r2_*k_r/(w1*Mn*((r1+r2_*k_r)^2+(x_s1+x_s2*k_x)^2))-ks;
end