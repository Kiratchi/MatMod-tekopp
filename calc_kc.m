function [R_calc_kc,Sc2,Sc] = calc_kc(T_top,p)


%T_film = (T_top + p.T_air)/2;
Pr_air = (my_air(p.T_air).*cp_air(p.T_air))./k_air(p.T_air); % Pr ok in range to use in schimdt 
D_AB = 2.634 * 1 / p.P_tot * (T_top/293)^(3/2);
%D_AB_P= D_AB*101325
%ny_water=my_water(T_top)/rho_water(T_top);
Sc = my_air(p.T_air) / (D_AB * rho_air(p.T_air));

% Sc1 = my_air(p.T_air) / (D_AB * rho_air(p.T_air))
Sc2 = my_air(T_top) / (D_AB * rho_air(T_top));

R_calc_kc = calc_h_l2air(T_top,p.r_inner*2)/(rho_air(p.T_air) * cp_air(p.T_air))*(Pr_air/Sc2)^(2/3);
end
%page 524 schimdt range 0.6 < Sc < 2500
%page 526 demonstrates the kc equation 
%page 393 eq. (22-41)  for D_AB and values used appendix table J.1 page 643

% there are multiple corrections compared to previous versions