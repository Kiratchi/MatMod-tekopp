
function q_values_calc(T_out_cup, T_top, p)
    
    global t_current time_values;
    time_values = [time_values, t_current];
    q_values = [q_rad_side(T_out_cup,p), q_rad_top(T_top,p), q_evap_top(T_top,p), q_top2air(T_top,p), q_glass2air(T_out_cup,p)];

    fileID = fopen('q_values.csv','a');
    if ~exist('q_values.csv', 'file')
        fprintf(fileID, '%s,%s,%s,%s,%s,%s\n', 'Time', 'q_rad_side', 'q_rad_top', 'q_evap_top', 'q_top2air', 'q_glass2air');
    end
    fprintf(fileID, '%f,%f,%f,%f,%f,%f\n', t_current, q_values);
    fclose(fileID); 

    relq_values=(q_values./sum(q_values));
    fileID = fopen('relq_values.csv','a');
    if ~exist('relq_values.csv', 'file')
        fprintf(fileID, '%s,%s,%s,%s,%s,%s\n', 'Time', 'relq_rad_side', 'relq_rad_top', 'relq_evap_top', 'relq_top2air', 'relq_glass2air');
    end
    fprintf(fileID, '%f,%f,%f,%f,%f,%f\n', t_current, relq_values);
    fclose(fileID);
end





