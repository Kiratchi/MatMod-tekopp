function final_q_glass = q_glass(T_out_cup,T_in_cup,p)

R_glass =  p.thickness_glass/p.A_side_ln*p.k_glass;
final_q_glass = (T_out_cup - T_in_cup)/R_glass;

end

