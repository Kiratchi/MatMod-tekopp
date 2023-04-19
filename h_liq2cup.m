function Resultat = h_liq2cup(T_in_cup, T_out_cup, C_L)
    % C_L is the characteristic length which in this case will be height to
    % which the water reaches in the container
    T_film = (T_in_cup + T_utsida)/2;
    konst_A=1.522; konst_B=243; my_water = konst_A*exp(konst_B/T_film);
    rho_water = 1000; % kg/m^3
    ny_water = my_water/rho_water;
    k_water = 0.598; % W / m K vid 20'C
    
    Pr_water =  50000/(T_film^2 + 155* T_film + 3700); % källa https://www.tec-science.com/mechanics/gases-and-liquids/prandtl-number/
    
    Ra_L = Pr_water * g * (1 / T_Abs) * (T_Abs - T_out_cup) * C_L^3 / (ny_water^2);
    
    Nu_L = 0.68 + 0.663*Ra_L^1/4/(1+(0.492/Pr_water)^(9/16))^(4/9) % om Ra_L <= 10^8
    
    Resultat = (k_water / C_L) * Nu_L
end 


