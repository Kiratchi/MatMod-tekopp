function n_A = calc_n_A(T_top,p)
C_air = p.RH * p_water(p.T_air)/(p.R*p.T_air);
C_surf = p_water(T_top)/(p.R*T_top);
n_A = calc_kc(Top,p)*p.A_top_l * p.M_water * (C_surf - C_air);
C_tot_surf = p.p_tot / (p.R*T_top);
mol_frac = C_surf / C_tot_surf;
n_A = n_A/(1-mol_frac);
end