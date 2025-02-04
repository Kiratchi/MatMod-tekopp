function q_comparer_model3(p,comp)
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
legend('q-rad', 'q-evap-top', 'q-top2air', 'q-glass2air');



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

    dTMdt(1) = -1/(cp_water(T_l)*rho_water(T_l)*p.volume_l)*(q_rad_side(T_l,p) + q_rad_top(T_l,p) + q_evap_top(T_l, p,0.5) + q_top2air(T_l,p,0.5) + q_glass2air(T_l,p,0.5) );
    dTMdt(2) = -calc_n_A(T_l, p,0.5);
    if (comp==1)
       
        time_values = [time_values, t_current];
        q_values_calc(T_l, T_l, p);
    end
end
end


function q_values_calc(T_out_cup, T_top, p)
    
    global t_current time_values;
    time_values = [time_values, t_current];
    q_values = [q_rad_side(T_out_cup,p) + q_rad_top(T_top,p), q_evap_top(T_top,p,0.5), q_top2air(T_top,p,0.5), q_glass2air(T_out_cup,p,0.5)];

    fileID = fopen('q_values.csv','a');
    if ~exist('q_values.csv', 'file')
        fprintf(fileID, '%s,%s,%s,%s,%s\n', 'Time', 'q_rad', 'q_evap_top', 'q_top2air', 'q_glass2air');
    end
    fprintf(fileID, '%f,%f,%f,%f,%f\n', t_current, q_values);
    fclose(fileID); 

    relq_values=(q_values./sum(q_values));
    fileID = fopen('relq_values.csv','a');
    if ~exist('relq_values.csv', 'file')
        fprintf(fileID, '%s,%s,%s,%s,%s\n', 'Time', 'relq_rad', 'relq_evap_top', 'relq_top2air', 'relq_glass2air');
    end
    fprintf(fileID, '%f,%f,%f,%f,%f\n', t_current, relq_values);
    fclose(fileID);
end
