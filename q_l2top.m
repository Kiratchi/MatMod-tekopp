function final_q_l2top = q_l2top(T_l,T_top,p)
    C_innerdiam = 2*p.r_inner; % Charateristic diameter 
    h_l2top = calc_h_l2top(T_top,T_l,C_innerdiam);
    R_l2top = (p.A_side_l*h_l2top)^-1;
    final_q_l2top = (T_l - T_top)/R_l2top;
end

