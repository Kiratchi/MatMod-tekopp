function R_calc_kc = calc_kc(T_top,p)

T_film = (T_top + p.T_air)/2;
D_AB = 2.634 * 1 / 101325 * (T_film/293)^(3/2);
ny_water=my_water(T_film)/rho_water(T_film);
Sc = ny_water / D_AB;
R_calc_kc = calc_h_l2air(T_top,p.r_inner)/(rho_air(T_film) * cp_air(T_film))*(pr_water(T_film)/Sc)^(2/3);
end