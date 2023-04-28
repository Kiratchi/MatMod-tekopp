function R_h_l2top = calc_h_l2top(T_top,T_l,C_innerdiam) 
    T_film = (T_top + T_l)/2;
    g = 9.82; % m / s^2
    %konst_A=1.522; konst_B=243; my_water = konst_A*exp(konst_B/T_film); %undersök enhet i andrade equation
    
    % my_water = 2.414 * 10^(-5) * 10^(247.8 / (T_film - 140));  % "Andrade equation." The original equation was proposed by E.N. da C. Andrade in a 1934 publication, "A Theory of the Viscosity of Liquids.
    
    A=1.856*10^-11; %mPa s 
    B=4209; %K 
    C=0.04527; %1/K 
    D=-3.376*10^-5; %1/K^2 
    my_water = A*exp(B/T_film + C*T_film +D*T_film^2);


    rho_water = 1000; % kg/m^3
    %k_water = 0.598; % W / m K vid 20'C
    ny_water = my_water / rho_water;
    %Pr_water = 50000/(T_film^2 + 155* T_film + 3700);
    Cp_water = 4180; % J / kg K
    Pr_water = (my_water*Cp_water(T_film))/k_water(T_film);

    % beta_water = 0.000207; % 1/K vid 20'C kan behöva använda istället för 1/T_l nedan

    Ra_L = g * Pr_water * 1/T_l * (T_l - T_top) * (C_innerdiam/4)^3 / (ny_water^2); % vanligtvis C_L^3 men C_L = A/O= pi r^2/(pi d) => d/4 ty (d/4)^3
     if (Ra_L >= 10^7) && (Ra_L <= 10^11)
        Nu_L = 0.15 * Ra_L^(1/3);
        R_h_l2top = (k_water(T_film) / (C_innerdiam/4))*Nu_L;
     elseif (Ra_L >= 10^4) && (Ra_L <= 10^7)
         Nu_L = 0.54 * Ra_L^(1/4);
         R_h_l2top = (k_water(T_film) / (C_innerdiam/4)) * Nu_L;
     else
        %warning('korrelationen håller ej för då Ra_L är för stor eller liten')
        Ra_L = 10^6;
        Nu_L = 0.54 * Ra_L^(1/4);
        R_h_l2top = (k_water(T_film) / (C_innerdiam/4))*Nu_L;
     end
end 

