%% Assumtions
% Lumped liquid box with air film around
% Assumes all water and glass have the same temperature.
% Optimize film temperature
clc, clear
for i = 1:5
    clf(figure(i))
end
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

% Physical properties
p.M_water = 18.01528 /1000; % kg/mol
p.R = 8.314; %kg⋅m2⋅s−2⋅K−1⋅mol−1
p.RH = 0.22;
p.P_tot = 101325; %Pa

% Physical dimensions
p.r_inner = 7 *10^(-2); %m
p.thickness_glass = 0.005; %m
p.r_outer = p.r_inner + p.thickness_glass; %m
p.height = 9.5*10^(-2); %m
p.A_side_l= 2*pi*p.r_inner*p.height; %m^2
p.A_side_ln = 2*pi*p.height/log(p.r_inner/p.r_outer); %m^2
p.A_side_glass = 2*pi*p.r_outer*p.height; %m^2
p.A_top_l = pi*p.r_inner^2; %m^2
p.volume_l = 250 *10^(-6); %m^3 %


% Teperature properties
p.T_air = 273.15+20.6; %K

% Radiation properies
p.sftboltz_const = 6.676*10^-8; %W/m^2 K^4 (Fundamentals of Momentum, Heat and Mass Transfer)
p.emissitivity_glass = 0.94;
p.emissitivity_l = 0.97 ;   %https://journals.ametsoc.org/view/journals/apme/11/8/1520-0450_1972_011_1391_ldowse_2_0_co_2.xml
p.rad_glass_const= p.A_top_l*p.emissitivity_glass*p.sftboltz_const;
p.rad_l_const= p.A_top_l*p.emissitivity_l*p.sftboltz_const;

% Transfer coefficents
p.k_glass = 0.9; %J/smK


plot_time_solution(p, 273.15+80, 150/1000, [0 2500])
plot_small_data()


figure(1)
axis([0 2500 15 100])
title("Change of temperature")
xlabel("Time (s)")
ylabel("T (C)")
legend('Our solution','Exp 1','Exp 2', 'Exp 3', 'Room temp')

figure(2)
title("Change of mass")
xlabel("Time (s)")
ylabel("Mass (g)")
legend('Our solution','Exp 1','Exp 2', 'Exp 3')

figure(3)
q_comparer_model3(p,1)

% For ful model
% saveas(figure(1),[pwd '/figures/T_over_t_m2'],'png')
% saveas(figure(2),[pwd '/figures/m_over_t_m2'],'png')
% saveas(figure(3),[pwd '/figures/q_compare_m2'],'png')


% For reduced model
saveas(figure(1),[pwd '/figures/T_over_t_m3'],'png')
saveas(figure(2),[pwd '/figures/m_over_t_m3'],'png')
saveas(figure(3),[pwd '/figures/q_compare_m3'],'png')


function dTMdt = derivate(p,TM_l,film_frac_side, film_frac_top)
    T_l = TM_l(1);
    T_M = TM_l(2);

    %For full model
    %dTMdt(1) = -1/(cp_water(T_l)*rho_water(T_l)*p.volume_l)*(q_rad_side(T_l,p) + q_rad_top(T_l,p) + q_evap_top(T_l, p,film_frac_top) + q_top2air(T_l,p,film_frac_top) + q_glass2air(T_l,p,film_frac_side) )
    
    %For reduced model
    dTMdt(1) = -1/(cp_water(T_l)*rho_water(T_l)*p.volume_l)*(q_evap_top(T_l, p,film_frac_top) + q_glass2air(T_l,p,film_frac_side) );
    

    dTMdt(2) = -calc_n_A(T_l, p,0.5);
end

function plot_time_solution(p, T_t0_l, M_t0_t, t_span) 
    f = @(t,TM) derivate(p,TM,0.5,0.5)';
    [t,y] = ode45(f, t_span, [T_t0_l M_t0_t]);
    T = y(:,1)-273.15;
    m = y(:,2)*1000

    figure(1)
    hold on
    plot(t, T,"black",'LineWidth',3)

    figure(2)
     hold on
    plot(t,m,"black",'LineWidth',3)
end

function sensitivity_analysis(p, T_t0_l, M_t0_t, t_span, span_side, span_top) 

    f = @(t,TM) derivate(p,TM,0.5,0.5)';
    [t,y] = ode45(f, t_span, [T_t0_l M_t0_t]);
    T = y(:,1)-273.15;
    m = y(:,2)*1000

    figure(4)
    hold on
    plot(t, T,"black",'LineWidth',3)

    figure(5)
     hold on
    plot(t,m,"black",'LineWidth',3)
end
