%% Assumtions
% Lumped liquid box with air film around
% Assumes all water and glass have the same temperature.
clc, clear, clf
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


derivate(p,[273.15+80, 100]);


% subplot(2,2,[1 2])
%plot_time_solution(p, 273.15+80, 125.81/1000, [0 158])
plot_time_solution(p, 273.15+80, 125.81/1000, [0 2500])
subplot(2,1,1)
plot_small_data()
legend('Our solution','Exp 1','Exp 2', 'Exp 3')

figure
q_comparer_model2(p,1)


function dTMdt = derivate(p,TM_l)
    T_l = TM_l(1);
    T_M = TM_l(2);
    %dTMdt(1) = -1/(cp_water(T_l)*rho_water(T_l)*p.volume_l)*(q_rad_side(T_l,p) + q_rad_top(T_l,p) + q_evap_top(T_l, p) + q_top2air(T_l,p) + q_glass2air(T_l,p) )
    dTMdt(1) = -1/(cp_water(T_l)*rho_water(T_l)*p.volume_l)*(q_evap_top(T_l, p) + q_glass2air(T_l,p) )
    

    dTMdt(2) = -calc_n_A(T_l, p);
end


function plot_time_solution(p, T_t0_l, M_t0_t, t_span) 
    f = @(t,TM) derivate(p,TM)';
    [t,y] = ode45(f, t_span, [T_t0_l M_t0_t]);
    T = y(:,1)-273.15;
    m = y(:,2)*1000;
    subplot(2,1,1)
    hold on
    plot(t, T,'LineWidth',1.5)
    %plot(t, 0*t+p.T_air-273.15, '--', 'Color', '#808080')
    axis([t_span 0 100])
    title("Change of temperature")
    xlabel("Time (s)")
    ylabel("T (C)")
    %ylim([60 80])
    subplot(2,1,2)
    plot(t,m,'LineWidth',1.5)
    axis([t_span 0 M_t0_t*1.1*1000])
    title("Change of mass")
    xlabel("Time (s)")
    ylabel("Mass (g)")
%     ylim([120 130])
%     xlim([0 20])

end

