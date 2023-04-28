function R_hcup2air = calc_h_cup2air(T_out_cup, C_L) 
    g=9.82; % m/s^2
    T_air = 293.15; % air is assumed to be constant
    T_film = (T_air+T_out_cup)/2; 
    Pr_air = 10^9/(1.1*(T_air-273.15)^3-1200*(T_air-273.15)^2+322000*(T_air - 273.15)+1.393*10^9);%0.7309 vid 20'C källa https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
    
   
    %k_air = 25.8*10^-3 ;% W / m K vid 20'C
    ny_air = my_air(T_film) / rho_air(T_film);

    Ra_L = Pr_air * g * 1/T_air * (T_out_cup-T_air) * C_L^3 / (ny_air^2);
    if (Ra_L <= 10^8)

        Nu_L = 0.68 + 0.67 * (Ra_L^0.25)/(1 + (0.492/Pr_air)^(9/16))^(4/9);
    else
        warning("Ra_L utanför correlationsintervall")
        Nu_L = 0.68 + 0.663*Ra_L^1/4/(1+(0.492/pr_water(T_film))^(9/16))^(4/9);
    end 

    R_hcup2air = (k_air(T_film)/ C_L) * Nu_L;

end

