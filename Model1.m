%% Assumtions
% Lumped liquid box with air film around. Only conduction i film
% No radiation
clc, clear, clf

T_t0_l = 365.15+100 %K
T_inf_a = 365.15+20.6 %K
C_p_l = 4.18   %Heat capacity water J/kg/K   %Funktion av temperaturen
density_l = 997 %kg/m^3
v_l = 250 *10^(-6) %m^3 
a_l= 2*pi*(7 *10^(-2))*9.5*10^(-2) %m^2
k_air= 0.6  %J/smK 
delta = 5*10^(-3) % m^2 Heat film thicknes

tspan = [0 2];


[t,y] = ode45(@(t,T) a_l*k_air/(C_p_l*density_l*v_l*delta)*(T_inf_a-T), tspan, T_t0_l);
plot(t,y)