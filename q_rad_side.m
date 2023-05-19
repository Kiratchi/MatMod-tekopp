function final_q_rad_side = q_rad_side(T_out_cup,p)
    rad_glass_const= p.A_top_l*p.emissitivity_glass*p.sftboltz_const;
    final_q_rad_side = rad_glass_const*(T_out_cup^4-p.T_air^4);

end

