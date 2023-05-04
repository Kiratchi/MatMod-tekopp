function R_calc_kc = calc_kc(T_top,p)


%T_film = (T_top + p.T_air)/2;
p.p_tot
D_AB = 2.634 * 1 / p.p_tot * (p.T_air/293)^(3/2)
%ny_water=my_water(T_top)/rho_water(T_top);
Sc = my_water(T_top) / (D_AB * rho_water(T_top))
R_calc_kc = calc_h_l2air(T_top,p.r_inner*2)/(rho_air(p.T_air) * cp_air(p.T_air))*(pr_water(p.T_air)/Sc)^(2/3);
end
%page 524 schimdt range 0.6 < Sc < 2500


