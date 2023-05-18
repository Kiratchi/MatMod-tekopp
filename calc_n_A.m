function n_A = calc_n_A(T_top,p,film_frac_top)
%T_film= (T_top+p.T_air)/2;
T_film = p.T_air*film_frac_top + T_top*(1-film_frac_top);

C_air = p.RH * p_water(p.T_air)/(p.R*p.T_air);
C_surf = p_water(T_film)/(p.R*T_film);
n_A = calc_kc(T_film,p,film_frac_top)*p.A_top_l * p.M_water * (C_surf - C_air);
C_tot_surf = p.P_tot / (p.R*T_film);
mol_frac = C_surf / C_tot_surf;
n_A = n_A/(1-mol_frac);
end