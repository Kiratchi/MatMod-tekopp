clc, clear,clf
 warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

data_large_beaker_1 = read_data('large_beaker_1.txt');
data_large_beaker_2 = read_data('large_beaker_2.txt');
data_large_beaker_3 = read_data('large_beaker_3.txt');
data_medium_beaker_1 = read_data('medium_beaker_1.txt');
data_medium_beaker_2 = read_data('medium_beaker_2.txt');
data_medium_surface_1 = read_data('medium_surface_1.txt');
data_medium_surface_2 = read_data('medium_surface_2.txt');


% hold on
% plot(data_large_beaker_1.T, data_large_beaker_1.t)
% plot(data_large_beaker_2.T, data_large_beaker_2.t)
% plot(data_large_beaker_3.T, data_large_beaker_3.t)
% 
% xlabel("t (s)")
% ylabel("T (C)")
% legend("exp. 1","exp. 2","exp. 3")
% title("Large beaker")


plot_small_data()



