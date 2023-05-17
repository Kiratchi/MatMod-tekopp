function final_q_rad_side = q_rad_side(T_out_cup,p)
    
    final_q_rad_side = p.rad_glass_const*(T_out_cup^4-p.T_air^4);

end

