function final_q_evap_top = q_evap_top(T_top,p,film_frac_top)
final_q_evap_top = calc_n_A(T_top,p,film_frac_top)*dHvap_water(T_top);
end

