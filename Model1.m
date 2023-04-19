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
p.emissitivity_l = 0.97    %https://journals.ametsoc.org/view/journals/apme/11/8/1520-0450_1972_011_1391_ldowse_2_0_co_2.xml


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

tspan = [0 100];
T_t0_l = 365.15+100; %K


[t,y] = ode45(@(t,T) derivate(p,T), tspan, T_t0_l);
plot(t,y)

function h = heat_transfer(p, function_flag)
    
    % function_flag: 1 för enbart toppen av cylindern, 2 för sidan av cylindern
    
    if function_flag == 1 
        h = 0.6 * (9.82 * 2.1 * 10^-4 * (p.T_area - p.T_air))^(1/4)*p.A_top_l/4;
    elseif function_flag == 2
        h = 0.6 * (9.82 * 2.1 * 10^-4 * (p.T_area - p.T_air))^(1/4)*p.A_top_l^2/4*p.height;
    else
        error('Ogiltlig function_flag: välj antingen 1 eller 2');
    end
end

function dTdt = derivate(p,T_l)
    h_l2air = 0.6; %Fixa utryck
    h_cup2air = 0.6; %Fixa utryck
    
    h_l2cup = 0; %Fixa utryck
    h_l2surf = 0; %Fixa utryck

    %R_l2glass = (p.A_side_l*h_l2glass)^-1;
    R_l2glass = 0 %Tilfäligt tills vi har h_l2glass
    R_glass =  p.thickness_glass/p.A_side_ln*p.k_glass;
    R_glass2air = (p.A_side_glass*h_cup2air)^-1;
    
    %R_l2surf = (p.A_side_l*h_l2surf)^-1;
    R_l2surf = 0 %Tillfälligt tills vi har h_l2surf
    R_surf2air = (p.A_top_l*h_l2air)^-1;


    T_surf_cup = (T_l*R_glass^-1 + p.T_air*R_glass2air^-1) / (R_glass^-1 + R_glass2air^-1); 
    T_surf_l = T_l;

    dQdt_side = (T_l-p.T_air)/(R_l2glass + R_glass + R_glass2air)
    dQdt_top = (T_l-p.T_air)/(R_l2surf + R_surf2air)

    dQdt_rad_side = p.A_side_glass*p.emissitivity_l*p.sftboltz_const*(T_surf_cup^4-p.T_air^4);  
    dQdt_rad_top = p.A_top_l*p.emissitivity_glass*p.sftboltz_const*(T_surf_l^4-p.T_air^4);

    dTdt = -1/(p.C_p_l*p.density_l*p.volume_l)*(dQdt_top + dQdt_side + dQdt_rad_top + dQdt_rad_side);
    
end
