function memristorplot(Data, save)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for x = 1:numel(Data)
        savefileroot = convertStringsToChars(Data(x).Filepath);
        savefileroot = savefileroot(1:end-4);
    
        tosave7 = figure;
        ax = gca;
        hold on
        plot(Data(x).Vout, Data(x).I, '-k');
        xlabel('Voltage (V)');
        ylabel('Current (A)');
       % axis([-4.1 2.6 -8E-3 8E-3]);
        yl = ylim(ax);
        axis(ax, 'tight')
        ylim(ax, yl)
        xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.05)
        box on
        set(gcf, 'Position', [100, 100, 600, 450])
        set(gca, 'FontSize', 18)
        hold off
        if save == 1
            saveas(tosave7,strcat(savefileroot,'_current.png'));
            saveas(tosave7,strcat(savefileroot,'_current.fig'));
        end
        
        tosave8 = figure;
        ax = gca;
        hold on
        plot(Data(x).Vout, Data(x).Iabs, '-k');
        xlabel('Voltage (V)');
        ylabel('Current (A)');
%         axis([-4.1 2.6 1.1E-8 1.4E-2]);
        yl = ylim(ax);
        axis(ax, 'tight')
        ylim(ax, yl)
        xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.05)
        box on
        set(gca, 'YScale', 'log');
        set(gcf, 'Position', [200, 100, 600, 450])
        set(gca, 'FontSize', 18)
        hold off
        if save == 1
            saveas(tosave8,strcat(savefileroot,'_LogCurrent.png'));
            saveas(tosave8,strcat(savefileroot,'_LogCurrent.fig'));
        end

        tosave9 = figure;
        ax = gca;
        hold on
        plot(Data(x).Vout, Data(x).R, '-k');
        xlabel('Voltage (V)');
        ylabel('Resistance (\Omega)');
        %axis([-4 2.5 2E2 inf]);
        yl = ylim(ax);
        axis(ax, 'tight')
        ylim(ax, yl)
        xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.05)
        box on
        set(gca, 'YScale', 'log');
        set(gcf, 'Position', [300, 100, 600, 450])
        set(gca, 'FontSize', 18)
        hold off
        if save == 1
            saveas(tosave9,strcat(savefileroot,'_Resistance.png'));
            saveas(tosave9,strcat(savefileroot,'_Resistance.fig'));
        end
        
        tosave9 = figure;
        ax = gca;
        hold on
        plot(Data(x).Vout, Data(x).R, '-k');
        xlabel('Voltage (V)');
        ylabel('Resistance (\Omega)');
        %axis([-4 2.5 2E2 inf]);
        yl = ylim(ax);
        axis(ax, 'tight')
        ylim(ax, yl)
        xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.05)
        box on
%         set(gca, 'YScale', 'log');
        set(gcf, 'Position', [300, 100, 600, 450])
        set(gca, 'FontSize', 18)
        hold off
        if save == 1
            saveas(tosave9,strcat(savefileroot,'_ResistanceLinear.png'));
            saveas(tosave9,strcat(savefileroot,'_ResistanceLinear.fig'));
        end
      
end
close all
end

