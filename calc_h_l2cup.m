function R_h_l2cup = calc_h_l2cup(T_in_cup, T_out_cup, C_L)
    % C_L is the characteristic length which in this case will be height to
    % which the water reaches in the container
    g = 9.82; %m/s^2 
    %  T_Abs = 273.15; % K
    T_film = (T_in_cup + T_out_cup)/2;
    %konst_A=1.522; konst_B=243; my_water = konst_A*exp(konst_B/T_film);

    A=1.856*10^-11/1000; %mPa s 
    B=4209; %K 
    C=0.04527; %1/K 
    D=-3.376*10^-5; %1/K^2 
    



    my_water = A*exp(B./T_film + C*T_film +D*T_film.^2); % "Andrade equation." The original equation was proposed by E.N. da C. Andrade in a 1934 publication, "A Theory of the Viscosity of Liquids.

    %rho_water = 1000; % kg/m^3
    ny_water = my_water/rho_water(T_film);
   % W / m K vid 20'C
    
    
    Pr_water = my_water * cp_water(T_film)/k_water(T_film);
    
    Ra_L = Pr_water * g * 1/T_film * (T_out_cup - T_in_cup) * C_L^3 / (ny_water^2);
    
    Nu_L = 0.68 + 0.663*Ra_L^1/4/(1+(0.492/Pr_water)^(9/16))^(4/9); % om Ra_L <= 10^8
    
    R_h_l2cup = (k_water / C_L) * Nu_L;
end 


