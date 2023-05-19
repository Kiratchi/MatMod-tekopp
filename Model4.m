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


film_frac_top = parameter_estimation(p, 1)
plot_time_solution(p, 273.15+80, 150/1000, [0 2500], film_frac_top)
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

% For best optim
saveas(figure(1),[pwd '/figures/T_over_t_optimized'],'png')
saveas(figure(2),[pwd '/figures/m_over_t_optimized'],'png')

% For T optim
% saveas(figure(1),[pwd '/figures/T_over_t_optimT'],'png')
% saveas(figure(2),[pwd '/figures/m_over_t_optimT'],'png')

% For m optim
% saveas(figure(1),[pwd '/figures/T_over_t_optimm'],'png')
% saveas(figure(2),[pwd '/figures/m_over_t_optimm'],'png')

% For side
% sensitivity_analysis(p, 273.15+80, 150/1000, [0 2500], [0.1 0.9], "film frac side")
% saveas(figure(4),[pwd '/figures/sensitivity_T_side'],'png')
% saveas(figure(5),[pwd '/figures/sensitivity_m_side'],'png')


% For top
% clf(figure(4))
% clf(figure(5))
% sensitivity_analysis(p, 273.15+80, 150/1000, [0 2500], [0.1 0.9], "film frac top")
% saveas(figure(4),[pwd '/figures/sensitivity_T_top'],'png')
% saveas(figure(5),[pwd '/figures/sensitivity_m_top'],'png')




function dTMdt = derivate(p,TM_l,film_frac_side,film_frac_top)
    T_l = TM_l(1);
    T_M = TM_l(2);

    %For full model
    %dTMdt(1) = -1/(cp_water(T_l)*rho_water(T_l)*p.volume_l)*(q_rad_side(T_l,p) + q_rad_top(T_l,p) + q_evap_top(T_l, p,film_frac_top) + q_top2air(T_l,p,film_frac_top) + q_glass2air(T_l,p,film_frac_side) )
    
    %For reduced model
    dTMdt(1) = -1/(cp_water(T_l)*rho_water(T_l)*p.volume_l)*(q_rad_side(T_l,p) + q_rad_top(T_l,p) + q_evap_top(T_l, p,film_frac_top) + q_glass2air(T_l,p,film_frac_side) );
    

    dTMdt(2) = -calc_n_A(T_l, p,film_frac_top);
end

function plot_time_solution(p, T_t0_l, M_t0_t, t_span, film_frac_top) 
    f = @(t,TM) derivate(p,TM, 0.5, film_frac_top)';
    [t,y] = ode45(f, t_span, [T_t0_l M_t0_t]);
    T = y(:,1)-273.15;
    m = y(:,2)*1000;

    figure(1)
    hold on
    plot(t, T,"black",'LineWidth',3)

    figure(2)
     hold on
    plot(t,m,"black",'LineWidth',3)
end

function sensitivity_analysis(p, T_t0_l, M_t0_t, t_span, span, type) 
    disp(type)
    light = [247, 241, 59]/255;
    dark = [59, 4, 51]/255;
    gradient = @(frac) light*frac + dark*(1-frac);

    
    f = @(t,TM) derivate(p,TM,0.5, 0.5)';
    [t,y] = ode45(f, t_span, [T_t0_l M_t0_t]);
    T = y(:,1)-273.15;
    m = y(:,2)*1000;

    figure(4)
    hold on
    plot(t, T,'LineWidth',3,'color', 'black')

    figure(5)
    hold on
    plot(t,m,'LineWidth',3,'color', 'black')
    
    for film_frac = linspace(span(1), span(2),20)
        if type =="film frac top"
            f = @(t,TM) derivate(p,TM,0.5, film_frac)';
        elseif type == "film frac side"
            f = @(t,TM) derivate(p,TM, film_frac,0.5)';
        end
        [t,y] = ode45(f, t_span, [T_t0_l M_t0_t]);
        T = y(:,1)-273.15;
        m = y(:,2)*1000;

        figure(4)
        hold on
        plot(t, T,'LineWidth',0.7,'color',gradient(film_frac))

        figure(5)
        hold on
        plot(t,m,'LineWidth',0.7,'color',gradient(film_frac))
        if film_frac == span(1) | film_frac == span(2)
            disp("T =" +  num2str(T(end)) + "  M =" +  num2str(m(end)))
        end
    end
    figure(4)
    title("Effect on temperature of changing " + type)
    xlabel("Time (s)")
    ylabel("Temperature (C)")

    figure(5)
    title("Effect on mass of changing " + type)
    xlabel("Time (s)")
    ylabel("Mass (g)")
end

function cost = param_costfunc(p, film_frac_top, data)
    cost = 0;
    for i = 1:3
        f = @(t,TM) derivate(p,TM,0.5, film_frac_top)';
        sol = ode45(f, [0 2500], [data{i}.T(1)+273.15 data{i}.m(1)/1000]);
        for j = 1:length(data{i}.t)
            y =  deval(sol,data{i}.t(j));
            T = y(1);
            m = y(2)*1000;
            % Only T
%             cost = cost + (T-data{i}.T(j)-273.15)^2;

            %Only m
%             cost = cost + (m-data{i}.m(j))^2;
            
            % Both
            cost = cost + ((T-data{i}.T(j)-273.15 )/(max(data{i}.T)-min(data{i}.T)))^2;
            cost = cost + ((m-data{i}.m(j) )/(max(data{i}.m)-min(data{i}.m)))^2;

        end
    end
end

function par_est = parameter_estimation(p, par_span)
    data{1} = read_data('small_beaker_1.txt');
    data{2} = read_data('small_beaker_2.txt');
    data{3} = read_data('small_beaker_3.txt');

    f = @(x) param_costfunc(p,x,data);
    par_est = fminbnd(f,0.1,0.9);
end
