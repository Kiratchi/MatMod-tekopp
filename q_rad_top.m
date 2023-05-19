function final_q_rad_top = q_rad_top(T_top,p)
    rad_l_const= p.A_top_l*p.emissitivity_l*p.sftboltz_const;
    final_q_rad_top = rad_l_const*(T_top^4-p.T_air^4);
    
end

