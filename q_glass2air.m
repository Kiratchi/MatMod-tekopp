function final_q_glass2air = q_glass2air(T_out_cup,p)
    
    C_L = p.height; % karaktäristisk längd 
    h_cup2air = calc_h_cup2air(T_out_cup, C_L);
    R_glass2air = (p.A_side_glass*h_cup2air)^-1;
    final_q_glass2air = (T_out_cup - p.T_air)/R_glass2air;

end

