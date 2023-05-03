function plot_small_data(Includelabels)
    data_small_beaker_1 = read_data('small_beaker_1.txt');
    data_small_beaker_2 = read_data('small_beaker_2.txt');
    data_small_beaker_3 = read_data('small_beaker_3.txt');
    data_small_surface_1 = read_data('small_surface_1.txt');
    data_small_surface_2 = read_data('small_surface_2.txt');
    data_small_surface_3 = read_data('small_surface_3.txt');

    hold on
    plot(data_small_beaker_1.T, data_small_beaker_1.t,'LineWidth',1.5)
    %plot(data_small_surface_1.T, data_small_surface_1.t,'b')

    plot(data_small_beaker_2.T, data_small_beaker_2.t,'LineWidth',1.5)
    %plot(data_small_surface_2.T, data_small_surface_2.t,'r')

    plot(data_small_beaker_3.T, data_small_beaker_3.t,'LineWidth',1.5)
    %plot(data_small_surface_3.T, data_small_surface_3.t,'k')
    plot(linspace(0,2500), ones(100)*20.6, '--', 'Color', '#808080')
    try
    if Includelabels ~=  false
        axis([0 2500 0 100])
        xlabel("t (s)")
        ylabel("T (C)")
        legend("exp. 1 beaker","exp. 1 surface","exp. 2 beaker","exp. 2 surface","exp. 3 beaker","exp. 3 surface")
        title("Small")
    end
    end
    
end