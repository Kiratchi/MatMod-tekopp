function R_h_l2air = calc_h_l2air(T_top,C_innerdiam)
    T_air = 273.15+20.6; % K
    T_film = (T_top + T_air)/2;
    %Pr_air = 10^9/(1.1*T_air^3-1200*T_air^2+322000*T_air+1.393*10^9); %0.7309 vid 20'C källa https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
    g=9.82; % m/s^2
    ny_air = my_air(T_film) / rho_air(T_film);


    Pr_air = (my_air(T_air).*cp_air(T_air))./k_air(T_air);
    Ra_L = g * Pr_air *(T_film - T_air) *1/T_air * (C_innerdiam/4)^3 / ((ny_air^2));

    if (Ra_L >= 10^7) && (Ra_L <= 10^11)
        Nu_L = 0.15 * Ra_L^(1/3);
        R_h_l2air = (k_air(T_film) / (C_innerdiam/4))*Nu_L;
     elseif (Ra_L >= 10^4) && (Ra_L <= 10^7)
        Nu_L = 0.54 * Ra_L^(1/4);
        R_h_l2air = (k_air(T_film) / (C_innerdiam/4)) * Nu_L;
    elseif (290 <= Ra_L) && (Ra_L <= 3.3*10^5)
        Nu_L = 0.78+0.82*Ra_L^(1/5); % källa: Martorell et al (2003) 
        R_h_l2air = (k_air(T_film) / (C_innerdiam/4)) * Nu_L;
    else
        %warning('the correlation does not  hold for the interval')
        %Disregard the meaning on this part just a warning to make the code
        %run but sinc we interate it still has to work.
        Ra_L = 10^6;
        Nu_L = 0.54 * Ra_L^(1/4);
        R_h_l2air = (k_air(T_film) / (C_innerdiam/4))*Nu_L;
     end

end 

