function R_hcup2air = calc_h_cup2air(T_out_cup, C_L) 
    g=9.82; % m/s^2
    T_air = 293.15; % air is assumed to be constant
    T_film = (T_air+T_out_cup)/2; 
    Pr_air = 10^9/(1.1*T_air^3-1200*T_air^2+322000*T_air+1.393*10^9);%0.7309 vid 20'C källa https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
    
    my_air_calc = my_air(T_film); % kg / m s vid 20'C
    rho_air_calc = rho_air(T_film); %kg/m^3 vid 20'C
    %k_air = 25.8*10^-3 ;% W / m K vid 20'C
    ny_air = my_air_calc / rho_air_calc;

    Ra_L = Pr_air * g * 1/T_air * (T_out_cup-T_air) * C_L^3 / (ny_air^2);

    Nu_L = 0.68 + 0.67 * (Ra_L^0.25)/(1 + (0.492/Pr_air)^(9/16))^(4/9);

    R_hcup2air = (k_air(T_film)/ C_L) * Nu_L;

end

