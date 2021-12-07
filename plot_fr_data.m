function plot_fr_data(frequencies, r, theta,style,sigout_ch,legend_on,fig_num,plot_phase,grid_on)

if(sigout_ch == 2)
   colors = ['g' 'm']; 
elseif(sigout_ch == 1)
   colors = ['b' 'r'];
elseif(sigout_ch == 3)
    colors = ['b' 'g'];
elseif(sigout_ch == 4)
    colors = ['r' 'm'];
elseif(sigout_ch == 5)
    colors = ['k' 'c'];
elseif(sigout_ch == 6)
    colors = ['k' 'y'];
end
% Plot data
%clf
figure(fig_num);
if(plot_phase)
  subplot(2, 1, 1)
end
s = plot(frequencies, 20*log10(r(1, :)), style);
hold on;
set(s, 'LineWidth', 2)
set(s, 'Color', colors(1));
s = plot(frequencies, 20*log10(r(2, :)), style);
set(s, 'LineWidth', 2)
set(s, 'Color', colors(2));

if(grid_on == 1)
    grid on
end

set(gca,'fontweight','bold','linewidth',2,'fontsize',12);

if(legend_on)
    if( (sigout_ch ==1) || (sigout_ch ==2))
        legend('in axis', 'cross axis');
    elseif(sigout_ch == 3)
        legend('drv180/sns0', 'drv225/sns45');
    elseif(sigout_ch == 4)
        legend('drv180/sns45', 'drv225/sns0');
    end
end

xlabel('Frequency [Hz]')
ylabel('Amplitude [dBVrms]')

if(plot_phase)
    subplot(2, 1, 2)
    s = plot(frequencies, theta(1, :)*180/pi, style);
    hold on;
    set(s, 'LineWidth', 2)
    set(s, 'Color', colors(1));
    s = plot(frequencies, theta(2, :)*180/pi, style);
    set(s, 'LineWidth', 2)
    set(s, 'Color', colors(2));
    if(grid_on == 1)
        grid on
    end
    xlabel('Frequency [Hz]')
    ylabel('Phase [deg]')
    set(gca,'fontweight','bold','linewidth',2,'fontsize',12);
end

end