function plot_small_data()
    data_small_beaker_1 = read_data('small_beaker_1.txt');
    data_small_beaker_2 = read_data('small_beaker_2.txt');
    data_small_beaker_3 = read_data('small_beaker_3.txt');
    data_small_surface_1 = read_data('small_surface_1.txt');
    data_small_surface_2 = read_data('small_surface_2.txt');
    data_small_surface_3 = read_data('small_surface_3.txt');


    figure(1)
    hold on
    plot(data_small_beaker_1.t, data_small_beaker_1.T,'LineWidth',1.5)
    %plot(data_small_surface_1.T, data_small_surface_1.t,'b')

    plot(data_small_beaker_2.t, data_small_beaker_2.T,'LineWidth',1.5)
    %plot(data_small_surface_2.T, data_small_surface_2.t,'r')

    plot(data_small_beaker_3.t, data_small_beaker_3.T,'LineWidth',1.5)
    %plot(data_small_surface_3.T, data_small_surface_3.t,'k')
    plot(linspace(0,2500), ones(100)*20.6, '--', 'Color', '#808080')
    

    figure(2)
    hold on
    data_small_beaker_1.t(1)
    plot(data_small_beaker_1.t, data_small_beaker_1.m +150-data_small_beaker_1.m(1),'LineWidth',1.5)
    plot(data_small_beaker_2.t, data_small_beaker_2.m +150-data_small_beaker_2.m(1),'LineWidth',1.5)
    plot(data_small_beaker_3.t, data_small_beaker_3.m +150-data_small_beaker_3.m(1),'LineWidth',1.5)
end