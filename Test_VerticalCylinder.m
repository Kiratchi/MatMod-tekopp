T_air=293.15+20;
T_out_cup=293.15+80;
g=9.82;
T_film=(T_air+T_out_cup)/2;
D = 7 *10^(-2)*2;
C_L = 9.5*10^(-2);
Pr_air = 10^9/(1.1*(T_air-273.15)^3-1200*(T_air-273.15)^2+322000*(T_air - 273.15)+1.393*10^9);%0.7309 vid 20'C källa https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm

ny_air = my_air(T_film) / rho_air(T_film);

Ra_L = Pr_air * g * 1/T_air * (T_out_cup-T_air) * C_L^3 / (ny_air^2);
Gr_L= Ra_L / Pr_air;
if (D/C_L) >= (35/Gr_L^(1/4))
    disp('Correlation for vertical wall holds')
end


