function R_deltaC = calc_deltaC(T_l,p)
p_A = p_water(T_l) * p.RH;
y_A = p_A/p.p_tot;%molar fraction in gas phase 
x_A = y_A * p.p_tot/p_water(T_l); % molar fraction in liquid phase
R_deltaC = x_A - y_A;
end

