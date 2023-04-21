%% Assumtions
% Lumped liquid box with air film around
clc, clear, clf


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

% Teperature properties
p.T_air = 273.15+20.6; %K
p.C_p_l = 4.18;   %Heat capacity water J/kg/K 

% Radiation properies
p.sftboltz_const = 6.676*10^-8; %W/m^2 K^4 (Transportenboken)
p.emissitivity_glass = 0.94;
p.emissitivity_l = 0.97 ;   %https://journals.ametsoc.org/view/journals/apme/11/8/1520-0450_1972_011_1391_ldowse_2_0_co_2.xml
p.rad_glass_const= p.A_top_l*p.emissitivity_glass*p.sftboltz_const;
p.rad_l_const= p.A_top_l*p.emissitivity_l*p.sftboltz_const;

% Transfer coefficents
p.k_glass = 0.9; %J/smK


tspan = [0 4000];
T_t0_l = 273.15+100; %K


[t,y] = ode45(@(t,T) derivate(p,T), tspan, T_t0_l);
T = y-273.15;
plot(t,T)
xlabel("Time (s)")
ylabel(" Temp (C)")
T(end)


function dTdt = derivate(p,T_l)


    %T_out_cup = (T_l*R_glass^-1 + p.T_air*R_glass2air^-1) / (R_glass^-1 + R_glass2air^-1); 
    %T_surf_l = T_l;

    T_0s = [273.15+80, 273.15+80]; %Fixa bättre initialgissning
    min_side = @(x) costfunc_side_flow(p,T_l, x(1), x(2));
    x = fminsearch(min_side,T_0s);
    T_in_cup = x(1);
    T_out_cup= x(2);
  
    T_top_0 = 273.15+80; %Fixa bättre initialgissning
    A=[];
    b=[];
    min_top = @(x) costfunc_top_flow(p,T_l, x);
    T_top = fmincon(min_top,T_top_0,A,b,273.15+20.6,273.15+100);


    C_innerdiam = 2*p.r_inner; % karaktäristisk diameter 
    h_l2top = calc_h_l2top(T_top,T_l,C_innerdiam); %Fixa utryck
    R_l2top = (p.A_side_l*h_l2top)^-1;
    q_l2top = (T_l - T_top)/R_l2top;

    C_L = p.height; % karaktäristisk längd 
    h_l2cup = calc_h_l2cup(T_in_cup, T_out_cup,C_L); %Fixa utryck
    R_l2glass = (p.A_side_l*h_l2cup)^-1;
    q_l2glass = (T_l - T_in_cup)/R_l2glass;


    dTdt = -1/(p.C_p_l*p.density_l*p.volume_l)*(q_l2top + q_l2glass);
end

function f = costfunc_top_flow(p,T_l, T_top)
    C_innerdiam = 2*p.r_inner; % karaktäristisk diameter 
    h_l2air = calc_h_l2air(T_top,C_innerdiam); %Fixa utryck
    
    h_l2top = calc_h_l2top(T_top, T_l, C_innerdiam); %Fixa utryck

    R_l2top = (p.A_side_l*h_l2top)^-1;
    R_top2air = (p.A_top_l*h_l2air)^-1;
    
    q_l2top = (T_l - T_top)/R_l2top;
    q_top2air = (T_top - p.T_air)/R_top2air;
    q_rad_top = p.rad_l_const*(T_top^4-p.T_air^4);

    f = (q_l2top - q_top2air - q_rad_top)^2;
end

function f = costfunc_side_flow(p,T_l, T_in_cup, T_out_cup)
    C_L = p.height; % karaktäristisk längd 
    h_l2cup = calc_h_l2cup(T_in_cup, T_out_cup,C_L);

    h_cup2air = calc_h_cup2air(T_out_cup, C_L);
    
    R_l2glass = (p.A_side_l*h_l2cup)^-1;
    R_glass =  p.thickness_glass/p.A_side_ln*p.k_glass;
    R_glass2air = (p.A_side_glass*h_cup2air)^-1;

    q_l2glass = (T_l - T_in_cup)/R_l2glass;
    q_glass = (T_out_cup - T_in_cup)/R_glass;
    q_glass2air = (p.T_air - T_out_cup)/R_glass2air;
    q_rad_side = p.rad_glass_const*(T_out_cup^4-p.T_air^4);

    f = (q_rad_side + q_glass2air - q_glass)^2 + (q_l2glass - q_glass)^2;
end






