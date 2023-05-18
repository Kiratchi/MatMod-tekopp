function final_q_top2air= q_top2air(T_top,p,film_frac_top)

    C_innerdiam = 2*p.r_inner; % Charateristic diameter 
    h_l2air = calc_h_l2air(T_top,C_innerdiam,p,film_frac_top);
    R_top2air = (p.A_top_l*h_l2air)^-1;
    final_q_top2air = (T_top - p.T_air)/R_top2air;
    
end

