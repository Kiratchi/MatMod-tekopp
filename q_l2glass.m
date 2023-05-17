function final_q_l2glass = q_l2glass(T_l,T_in_cup,p)
    C_L = p.height; % Charecteristic length
    h_l2cup = calc_h_l2cup(T_in_cup, T_l,C_L);
    R_l2glass = (p.A_side_l*h_l2cup)^-1;
    final_q_l2glass = (T_l - T_in_cup)/R_l2glass;
end

