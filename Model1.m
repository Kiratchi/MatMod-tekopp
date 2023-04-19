%% Assumtions
% Lumped liquid box with air film around. Only conduction i film
% No radiation
clc, clear, clf

% Teperature properties
p.T_air = 365.15+20.6; %K
p.C_p_l = 4.18;   %Heat capacity water J/kg/K 

% Radiation properies
p.sftboltz_const = 6.676*10^-8 %W/m^2 K^4 (Transportenboken)
p.emissitivity_glass = 0.94
p.emissitivity_water = 0.97    %https://journals.ametsoc.org/view/journals/apme/11/8/1520-0450_1972_011_1391_ldowse_2_0_co_2.xml


% Transfer coefficents
p.k_glass = 0.9; %J/smK

% Physical dimensions
p.r_inner = 7 *10^(-2); %m
p.thickness_glass = 0.005; %m
p.r_outer = p.r_inner + p.thickness_glass; %m
p.height = 9.5*10^(-2); %m
p.A_side_l= 2*pi*p.r_inner*p.height; %m^2
p.A_side_ln = 2*pi*p.height/log(p.r_inner/p.r_outer); %m^2
p.A_side_glass = 2*pi*p.r_outer*p.height; %m^2
p.A_top_l = pi*p.r_inner^2; %m^2
p.volume_l = 250 *10^(-6); %m^3 

p.density_l = 997; %kg/m^3

tspan = [0 200];
T_t0_l = 365.15+100; %K


[t,y] = ode45(@(t,T) derivate(p,T), tspan, T_t0_l);
plot(t,y)

function dTdt = derivate(p,T_l)
    %h_l2air = 0.6*(9.82*beta * (T_surface-T_air))^(1/4)*
    %h_cup2air=0.6;
    h_l2air = 0.6;
    h_cup2air = 0.6;

    T_surf_cup = (p.A_side_ln*p.k_glass*T_l + p.A_side_glass*h_cup2air*p.T_air)...
        /(p.A_side_ln*p.k_glass + p.A_side_glass*h_cup2air);
    T_surf_l = T_l;

    dQdt_conv_side = p.A_side_ln*p.k_glass*(T_surf_cup-p.T_air); % Något fel här
    dQdt_conv_top = p.A_top_l*h_l2air*(T_surf_l-p.T_air);
    dQdt_conv_side = 0

    dQdt_rad_side = p.A_side_glass*p.emissitivity_water*p.sftboltz_const*(T_surf_cup^4-p.T_air^4);  
    dQdt_rad_top = p.A_top_l*p.emissitivity_glass*p.sftboltz_const*(T_surf_l^4-p.T_air^4);

    dTdt = -1/(p.C_p_l*p.density_l*p.volume_l)*(dQdt_conv_top + dQdt_conv_side + dQdt_rad_top + dQdt_rad_side);
    
end
