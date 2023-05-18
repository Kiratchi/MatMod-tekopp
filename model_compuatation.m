function model_compuatation(p,comp)
global time_values;
global t_current
time_values = [];
clearfile('relq_values.csv')
clearfile('q_values.csv')
T_t0_l=273.15+80;
M_t0_t=125.81/1000;
t_span=[0 2500];
f = @(t,TM) derivate(p,TM,1,t)';
t_current = t_span(1);
options = odeset('OutputFcn', @odeplot, 'OutputSel', 1,'RelTol',1e-4,'AbsTol',1e-4);
[t,y] = ode45(f, t_span, [T_t0_l M_t0_t], options);
T = y(:,1)-273.15;
m = y(:,2)*1000;
size(t);

data = csvread('relq_values.csv', 1, 0); 
size(data,1);
for i = 2:size(data, 2)
    plot(data(:, 1),data(:, i),'LineWidth',2);
    hold on;
end
xlabel('Time');
ylabel('Relative contribution of different q terms');
legend('q-rad-side', 'q-rad-top', 'q-evap-top', 'q-top2air', 'q-glass2air');



function status = odeplot(t,~,~)
  
  if ~isempty(t)  % Check if t is not empty
      t_current = t(end);
  end
  status = 0;  
end

function dTMdt = derivate(p, TM_l, comp,t)
   
    t_current = t; % Set the global variable to the current time

    T_l = TM_l(1);
    T_M = TM_l(2);
    [T_in_cup, T_out_cup] = t_finder_side(p, T_l);
    T_top = t_finder_top(p, T_l);

    dTMdt(1) = -1/(cp_water(T_l)*rho_water(T_l)*p.volume_l)*(q_l2top(T_l,T_top,p) + q_l2glass(T_l,T_in_cup,p));
    dTMdt(2) = -calc_n_A(T_top, p);
    if (comp==1)
       
        time_values = [time_values, t_current];
        q_values_calc(T_out_cup, T_top, p);
    end
end


%--------------------------------------------------------------------------
function [T_in_cup, T_out_cup] = t_finder_side(p,T_l)
    options = optimoptions('fmincon','Display', 'off');
    min_side = @(x) costfunc_side_flow(p,T_l, x(1), x(2));
    %disp("T_l =" + (T_l-273.15))
    [x,f_val] = fmincon(min_side,[T_l-5, T_l-4],[],[],[],[],[T_l-60 T_l-60],[T_l+20 T_l+20],[],options);
    T_in_cup = x(1);
    T_out_cup= x(2);
end

function T_top = t_finder_top(p,T_l)
    options = optimoptions('fmincon','Display', 'off');
    min_top = @(x) costfunc_top_flow(p,T_l, x);
    T_top = fmincon(min_top,T_l-5,[],[],[],[],273.15+20.6,273.15+100, [],options);
end

function f = costfunc_top_flow(p,T_l, T_top)

    f = (q_l2top(T_l,T_top,p) - q_top2air(T_top,p) - q_rad_top(T_top,p)-q_evap_top(T_top,p))^2;
end

function f = costfunc_side_flow(p,T_l, T_in_cup, T_out_cup)


    f = (q_rad_side(T_out_cup,p) + q_glass2air(T_out_cup,p) - q_glass(T_out_cup,T_in_cup,p))^2 + (q_l2glass(T_l,T_in_cup,p) - q_glass(T_out_cup,T_in_cup,p))^2;
end
end