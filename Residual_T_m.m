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


film_frac_top = parameter_estimation(p, 1);
plot_time_solution(p, 273.15+80, 150/1000, [0 2500], film_frac_top)
plot_small_data_res(p)




%---------------------------------------------------
% Residual plotts and calculatin R2 values 
%---------------------------------------------------
% Plot solution and get predicted data
[T_pred, m_pred, t] = plot_time_solution(p, 273.15+80, 150/1000, [0 2500], 0.5);

% Read experimental data from files
data1 = read_data('small_beaker_1.txt');
data2 = read_data('small_beaker_2.txt');
data3 = read_data('small_beaker_3.txt');

% Perform the operations for each data set
datasets = {data1, data2, data3};
R2_T = zeros(length(datasets), 1);
R2_m = zeros(length(datasets), 1);
% Initialize vectors for all experimental and predicted data
T_exp_all = [];
m_exp_all = [];
T_pred_all = [];
m_pred_all = [];

max_length = max([length(datasets{1}.T), length(datasets{2}.T), length(datasets{3}.T)]);

% Make all datasets the same length 
T_exp_padded = NaN(max_length, length(datasets));
m_exp_padded = NaN(max_length, length(datasets));


for i = 1:length(datasets)
    T_exp_padded(1:length(datasets{i}.T), i) = datasets{i}.T;
    m_exp_padded(1:length(datasets{i}.m), i) = datasets{i}.m;
end

% Calculate experimental error
exp_error_T = std(T_exp_padded, 'omitnan');
exp_error_m = std(m_exp_padded, 'omitnan');

for i = 1:length(datasets)
    data = datasets{i};

    % Get experimental data
    T_exp = data.T;
    m_exp = data.m;
    t_exp = data.t;

    % Interpolate the model predictions onto the same time 
    T_pred_interp = interp1(t, T_pred, t_exp);
    m_pred_interp = interp1(t, m_pred, t_exp);

    % Append to all data vectors
    T_exp_all = [T_exp_all; T_exp];
    m_exp_all = [m_exp_all; m_exp];
    T_pred_all = [T_pred_all; T_pred_interp];
    m_pred_all = [m_pred_all; m_pred_interp];
end

% Calculate R2 values for all data
SSE_T = sum((T_pred_all - T_exp_all).^2);
SST_T = sum((T_exp_all - mean(T_exp_all)).^2);
R2_T = 1 - SSE_T/SST_T;



SSE_m = sum((m_pred_all - m_exp_all).^2);
SST_m = sum((m_exp_all - mean(m_exp_all)).^2);
R2_m = 1 - SSE_m/SST_m;



% Print the R2 values and experimental errors
disp(['Temp R2: ' num2str(R2_T) ', Experimental error: ' num2str(exp_error_T)]);
disp(['Mass R2: ' num2str(R2_m) ', Experimental error: ' num2str(exp_error_m)]);

%---------------------------------------------------


function dTMdt = derivate(p,TM_l,film_frac_side,film_frac_top)
    T_l = TM_l(1);
    T_M = TM_l(2);

    %For full model
    %dTMdt(1) = -1/(cp_water(T_l)*rho_water(T_l)*p.volume_l)*(q_rad_side(T_l,p) + q_rad_top(T_l,p) + q_evap_top(T_l, p,film_frac_top) + q_top2air(T_l,p,film_frac_top) + q_glass2air(T_l,p,film_frac_side) )
    
    %For reduced model
    dTMdt(1) = -1/(cp_water(T_l)*rho_water(T_l)*p.volume_l)*(q_rad_side(T_l,p) + q_rad_top(T_l,p) + q_evap_top(T_l, p,film_frac_top) + q_glass2air(T_l,p,film_frac_side) );
    

    dTMdt(2) = -calc_n_A(T_l, p,film_frac_top);
end

function [T,m,t] = plot_time_solution(p, T_t0_l, M_t0_t, t_span, film_frac_top) 
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
            %cost = cost + (T-data{i}.T(j)-273.15)^2;

            %Only m
            %cost = cost + (m-data{i}.m(j))^2;
            
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

function plot_small_data_res(p)
    data_small_beaker_1 = read_data('small_beaker_1.txt');
    data_small_beaker_2 = read_data('small_beaker_2.txt');
    data_small_beaker_3 = read_data('small_beaker_3.txt');
    data_small_surface_1 = read_data('small_surface_1.txt');
    data_small_surface_2 = read_data('small_surface_2.txt');
    data_small_surface_3 = read_data('small_surface_3.txt');

    [T_pred, m_pred, t] = plot_time_solution(p, 273.15+80, 150/1000, [0 2500], 0.5);

    data_beakers = {data_small_beaker_1, data_small_beaker_2, data_small_beaker_3};

    figure(3)
    hold on

    for i = 1:3
        T_pred_interp = interp1(t, T_pred, data_beakers{i}.t);
        residuals_T = data_beakers{i}.T - T_pred_interp;  % Temperature residuals

        plot(data_beakers{i}.t, residuals_T, '*')
        title('Residuals of Temperature')
        xlabel('Time (s)')
        ylabel('Residuals')
    end

    legend('Experiment 1', 'Experiment 2', 'Experiment 3')

    yline(0,'--k')

    figure(4)
    hold on

    for i = 1:3
        m_pred_interp = interp1(t, m_pred, data_beakers{i}.t);
        residuals_m = (data_beakers{i}.m + 150 - data_beakers{i}.m(1)) - m_pred_interp;  % Mass residuals

        plot(data_beakers{i}.t, residuals_m, '*')
        title('Residuals of Mass')
        xlabel('Time (s)')
        ylabel('Residuals')
    end

    legend('Experiment 1', 'Experiment 2', 'Experiment 3')

    yline(0,'--k')
    
end



