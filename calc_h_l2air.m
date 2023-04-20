function R_h_l2air = calc_h_l2air(T_top,C_innerdiam)
    T_air = 293.15; % K
    T_film = (T_top + T_air)/2;
    Pr_air = 0.7309; % vid 20'C källa https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
    g=9.82; % m/s^2
    my_air = 1.825*10^-5; % kg / m s vid 20'C kan göras temperaturberoende 
    rho_air = 1.204; %kg/m^3 vid 20'C
    k_air = 25.8*10^-3; % W / m K vid 20'C
    ny_air = my_air / rho_air;

    Gr_L = g * (T_top - T_air) *1/T_air * (C_innerdiam/4)^3 / ((ny_air^2));
    Ra_L = Gr_L * Pr_air;
    Nu_L = 0.54 * Ra_L^(1/4); % implementera kanske en if för att kolla Ra_L och ändra korrelation.
    R_h_l2air = k_air/(C_innerdiam/4) * Nu_L;

end 

