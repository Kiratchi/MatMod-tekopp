%% Assumtions
% Lumped liquid box with air film around
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
q_comparer_model2(p,1)

figure(4)
plot_side_temp(p,p.T_air+5, 273.15+100)

figure(5)
plot_top_temp(p,p.T_air+5, 273.15+100)

% saveas(figure(1),[pwd '/figures/T_over_t_m1'],'png')
% saveas(figure(2),[pwd '/figures/m_over_t_m1'],'png')
% saveas(figure(3),[pwd '/figures/q_compare_m1'],'png')
% saveas(figure(4),[pwd '/figures/side_effect_m1'],'png')
% saveas(figure(5),[pwd '/figures/top_effect_m1'],'png')

function dTMdt = derivate(p,TM_l)
    T_l = TM_l(1);
    T_M = TM_l(2);
    [T_in_cup, T_out_cup] = t_finder_side(p, T_l);
    T_top = t_finder_top(p, T_l);    
  
    dTMdt(1) = -1/(cp_water(T_l)*rho_water(T_l)*p.volume_l)*(q_l2top(T_l,T_top,p) + q_l2glass(T_l,T_in_cup,p));
    dTMdt(2) = -calc_n_A(T_top, p,0.5);
end


function [T_in_cup, T_out_cup] = t_finder_side(p,T_l)
    options = optimoptions('fmincon','Display', 'off');
    min_side = @(x) (q_rad_side(x(2),p) + q_glass2air(x(2),p,0.5) - q_glass(x(2),x(1),p))^2 + (q_l2glass(T_l,x(1),p) - q_glass(x(2),x(1),p))^2;
    [x,f_val] = fmincon(min_side,[T_l-5, T_l-4],[],[],[],[],[T_l-60 T_l-60],[T_l+20 T_l+20],[],options);
    T_in_cup = x(1);
    T_out_cup= x(2);
end

function T_top = t_finder_top(p,T_l)
    options = optimoptions('fmincon','Display', 'off');
    min_top = @(x) (q_l2top(T_l,x,p) - q_top2air(x,p,0.5) - q_rad_top(x,p)-q_evap_top(x,p,0.5))^2;
    T_top = fmincon(min_top,T_l-5,[],[],[],[],273.15+20.6,273.15+100, [],options);
end

function plot_side_temp(p, T_low, T_high)
T_l_vector = linspace(T_low, T_high,100);
T_in_cup_vector = [];
T_out_cup_vector= [];
for T_l = T_l_vector
    [T_in_cup, T_out_cup] = t_finder_side(p, T_l);
    T_in_cup_vector = [T_in_cup_vector, T_in_cup];
    T_out_cup_vector = [T_out_cup_vector, T_out_cup];
end
disp("Avrage diffrence between T_in_cup and T_l is:" + sum(T_l_vector-T_in_cup_vector)/length(T_l_vector))
disp("Avrage diffrence between T_out_cup and T_l is:" + sum(T_l_vector-T_out_cup_vector)/length(T_l_vector))
disp("Avrage diffrence between T_out_cup and T_in_cup is:" + sum(T_in_cup_vector-T_out_cup_vector)/length(T_in_cup_vector))
plot(T_l_vector-273.15, T_in_cup_vector-273.15,'-','linewidth', 1.5)
hold on 
plot(T_l_vector-273.15, T_out_cup_vector-273.15,'--','LineWidth',1.5)
xlabel("T_l (C)")
ylabel("T_{cup} (C)")
legend("Outside", "Inside")
title("Correlation T_l & Side temp")
end


function plot_top_temp(p, T_low, T_high)
T_l_vector = linspace(T_low, T_high,100);
T_top_vector = [];
for T_l = T_l_vector
    T_top = t_finder_top(p, T_l);
    T_top_vector = [T_top_vector,T_top ];
end
disp("Avrage diffrence between T_top and T_l is: " + sum(T_l_vector-T_top_vector)/length(T_l_vector))
plot(T_l_vector-273.15, T_top_vector-273.15,'-','linewidth', 1.5)
hold on 
xlabel("T_l (C)")
ylabel("T_{top} (C)")
legend("Top")
title("Correlation T_l & top temp")
end

function plot_time_solution(p, T_t0_l, M_t0_t, t_span) 

    f = @(t,TM) derivate(p,TM)';
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

