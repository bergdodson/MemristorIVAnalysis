close all;
FigurePlot = 1; %tells whether the code will go ahead with programmed plots or just import data.
plotnum = 1; %determines which set of data will be plotted
quest1 = 'Do you want to select a new folder to import?';
answer1 = questdlg(quest1);
save = 1;
SetFilter = 150; %controls how much above background the instantaneous derivative needs to be to be considered a SET/RESET switch;
ResetFilter = 50; %controls how much above background the instantaneous derivative needs to be to be considered a SET/RESET switch;
Thresh = 1.0; %controls threshold for how program finds false switches that reverse themselves (this is a percentage of the 1st derivative).
switch  answer1
    case 'Yes'
        R0 = 0; %initializes initial resistance of device
        RatioMax = 0; %Reinitializes the max ratio found when reading all the data files
        Rmax = 0;%HRS for max ratio at 100 mV
        Rmin = 0; %LRS for max ratio at 100 mV
        RatioMax2 = 0; %Reinitializes the 2nd highest ratio found, otherwise previous program runs would influence this.
        Rmax2 = 0;
        Rmin2 = 0;
        path = uigetdir; %Brings up the GUI for selecting a folder to load everything from
        filepattern = fullfile(path, '*.txt'); %Sets the filepath and filetype that the directory search will return
        Filelist = dir(filepattern); %extracts all filenames in the given directory
        Filenames = cell(size(Filelist, 1), 1); %initializing 'Filenames' cell array
        Data = struct('Filename', [], 'Filepath', [], 'Vout', [], 'I', [], 'Iabs', [], 'R', [],'Vstart', [], 'Vend', [], 'SweepType', [], 'VsetRawIndex', [], 'VsetRaw', [], 'VsetDiff', [], 'VsetStatus', [], 'Vset', [], 'VresetRawIndex', [], 'VresetRaw', [], 'VresetDiff', [], 'VresetStatus', [], 'Vreset', [], 'R100F', [], 'R100R', [], 'R100ratio', [], 'Vmax', [], 'Vmin', []); %makes a 'struct' known as 'Data' that is the same dimensions of the # of files in the selected folder
        for c = 1:size(Filelist,1)
            Filenames(c) =  cellstr(Filelist(c).name); %makes a cell array out of extracted file names which are strings so the names can be sorted properly
        end

        Filenames = natsort(Filenames); %imported from matlab website, sorts filenames properly by the number so instead of 1, 10, 11, ... we get 1, 2, 3

        for c = 1:size(Filelist,1)
            Data(c).Filename = string(Filenames(c)); %changing sorted file name back to string
            Data(c).Filepath = strcat(path, '\', Data(c).Filename); %creating a full file path out of the strings
            A = importdata(Data(c).Filepath, '\t', 110); %Importing the data with tab delimiting and skipping the first 110 lines      
            Data(c).Vout = A.data(:,1);  %next few lines are just reading imported data into a struct array and processing it.
            Data(c).I = A.data(:,2);
            Data(c).Iabs = abs(A.data(:,2));
            Data(c).R = A.data(:,3);
            Data(c).Vstart = Data(c).Vout(1); %1st voltage of the list, so your starting sweep value
            Data(c).Vend = Data(c).Vout(size(Data(c).Vout, 1)); %reads last data point using last index to find ending point of voltage sweep
            if Data(c).Vstart == Data(c).Vend %If the software detects that a sweep is a full loop like starting and ending at the same point it will then split the data in half to analyze the forward and backwards sweeps
                Data(c).SweepType = "Double";
                a = size(Data(c).Vout, 1);%gives the total number of points
                b = a/2; % gives the half index value
                V1 = Data(c).Vout(1:b);%splitting the forward and reverse sweep into two different data sets for interpolation.
                R1 = Data(c).R(1:b);
                I1 = Data(c).I(1:b);
                D1 = diff(I1);%takes instantaneous derivatives of all the resistances to detect changes
                Avg1 = mean(rmoutliers(abs(D1))); %finds average value of forward sweep with spikes from switching removed
                [Junk, Data(c).VsetRawIndex] = findpeaks(abs(D1), 'MinPeakHeight', SetFilter*Avg1); %gives switching derivative values filtered into a junk variable 'Dump' and the indices they're found at.
                %This loop fills in the index and raw unfiltered SET voltage values detected
                for j = 1:size(Data(c).VsetRawIndex)
                    Data(c).VsetRaw(j) = V1(Data(c).VsetRawIndex(j));
                    Data(c).VsetDiff(j) = D1(Data(c).VsetRawIndex(j));
                end

                %The following loop checks all the detected switching values and filter them
                k = 0; % will be used to control new filtered Vset data set, and remove anomalies like back-switching
                for j = 1:size(Data(c).VsetRawIndex, 1)
                    if Data(c).VsetRaw(j) <=0                                                   %for checking if negative bias switch happened
                        Data(c).VsetStatus(j) = "Negative Bias switch";
                    elseif j == size(Data(c).VsetRawIndex, 1) && Data(c).VsetDiff(j) > 0        %for last switch in set if proper switch
                        k = k+1;
                        Data(c).VsetInd(k) = Data(c).VsetRawIndex(j);
                        Data(c).Vset(k) = Data(c).VsetRaw(j);
                        Data(c).VsetStatus(j) = "Normal";
                    elseif j == size(Data(c).VsetRawIndex, 1) && Data(c).VsetDiff(j) <= 0        %for last switch in set if negative
                        Data(c).VsetStatus(j) = "Negative";
                    elseif Data(c).VsetDiff(j) > 0 && Data(c).VsetDiff(j+1)>0                    %checking if a switch is positive and followed by another positive switch afterwards
                        k = k+1;
                        Data(c).VsetInd(k) = Data(c).VsetRawIndex(j);
                        Data(c).Vset(k) = Data(c).VsetRaw(j);
                        Data(c).VsetStatus(j) = "Normal";
                    elseif Data(c).VsetDiff(j) > 0 && Data(c).VsetDiff(j+1)<=0 && Data(c).VsetRaw(j+1)-Data(c).VsetRaw(j) >= 0.5 %checks if distance between pos/neg switch is sufficiently large to consider it not a flipflop
                        k = k+1;
                        Data(c).VsetInd(k) = Data(c).VsetRawIndex(j);
                        Data(c).Vset(k) = Data(c).VsetRaw(j);
                        Data(c).VsetStatus(j) = "Normal";
                    elseif Data(c).VsetDiff(j) > 0 && Data(c).VsetDiff(j+1)<=0                   %checking if a switch is positive but followed by a negative switch.
                        Perc = (Data(c).VsetDiff(j)+Data(c).VsetDiff(j+1))/Data(c).VsetDiff(j); %finds percentage difference between the two switches based on set threshold value 'Thresh'
                        if abs(Perc) <= Thresh                                                  % if under the threshold, the switch is noted as a reversed switch, but not added to the final filtered set
                            Data(c).VsetStatus(j) = "Net 0 flip";
                        elseif abs(Perc) > Thresh && Perc < 0
                            k = k+1;
                            Data(c).VsetInd(k) = Data(c).VsetRawIndex(j);
                            Data(c).Vset(k) = Data(c).VsetRaw(j);
                            Data(c).VsetStatus(j) = "Net RESET flip";
                        elseif abs(Perc) > Thresh && Perc > 0
                            k = k+1;
                            Data(c).VsetInd(k) = Data(c).VsetRawIndex(j);
                            Data(c).Vset(k) = Data(c).VsetRaw(j);
                            Data(c).VsetStatus(j) = "Net SET flip";
                        else
                            Data(c).VsetStatus(j) = "Error";
                        end
                    elseif Data(c).VsetDiff(j) < 0
                        Data(c).VsetStatus(j) = "Negative";
                    else
                        Data(c).VsetStatus(j) = "Error";
                    end
                end

                Data(c).R100F = max(interp1(V1, R1, 0.1, 'nearest'), interp1(V1, R1, -0.1, 'nearest')); %finding resistance at or near 0.100 V (positive or negative, it picks the biggest one)

                V2 = Data(c).Vout(b+1:a);
                R2 = Data(c).R(b+1:a);
                I2 = Data(c).I(b+1:a);
                D2 = diff(I2); 
                Avg2 = mean(rmoutliers(abs(D2))); %finds average value of forward sweep with spikes from switching removed
                [Junk, Data(c).VresetRawIndex] = findpeaks(abs(D2), 'MinPeakHeight', ResetFilter*Avg2); %gives switching voltages and the indices they're found at.
                for j = 1:size(Data(c).VresetRawIndex)
                    Data(c).VresetRaw(j) = V2(Data(c).VresetRawIndex(j));
                    Data(c).VresetDiff(j) = D2(Data(c).VresetRawIndex(j));
                end
                
                k = 0; % will be used to control new filtered Vset data set, and remove anomalies like back-switching
                for j = 1:size(Data(c).VresetRawIndex, 1)
                    if Data(c).VresetRaw(j) >=0
                        Data(c).VresetStatus(j) = "Positive Bias switch";
                    elseif j == size(Data(c).VresetRawIndex, 1) && Data(c).VresetDiff(j) > 0 %for last switch in set if proper switch
                        k = k+1;
                        Data(c).VresetInd(k) = Data(c).VresetRawIndex(j);
                        Data(c).Vreset(k) = Data(c).VresetRaw(j);
                        Data(c).VresetStatus(j) = "Normal";
                    elseif j == size(Data(c).VsetRawIndex, 1) && Data(c).VsetDiff(j) <= 0 %for last switch in set if negative
                        Data(c).VresetStatus(j) = "Negative";
                    elseif Data(c).VresetDiff(j) > 0 && Data(c).VresetDiff(j+1)>0  %checking if a switch is positive and followed by another positive switch afterwards
                        k = k+1;
                        Data(c).VresetInd(k) = Data(c).VresetRawIndex(j);
                        Data(c).Vreset(k) = Data(c).VresetRaw(j);
                        Data(c).VresetStatus(j) = "Normal";
                    elseif Data(c).VresetDiff(j) > 0 && Data(c).VresetDiff(j+1)<=0 && Data(c).VresetRaw(j+1)-Data(c).VresetRaw(j) <= -0.5
                        k = k+1;
                        Data(c).VresetInd(k) = Data(c).VresetRawIndex(j);
                        Data(c).Vreset(k) = Data(c).VresetRaw(j);
                        Data(c).VresetStatus(j) = "Normal";
                    elseif Data(c).VresetDiff(j) > 0 && Data(c).VresetDiff(j+1)<=0 %checking if a switch is positive but followed by a negative switch.
                        Perc = (Data(c).VresetDiff(j)+Data(c).VresetDiff(j+1))/Data(c).VresetDiff(j); %finds percentage difference between the two switches
                        if abs(Perc) <= Thresh % if under the threshold, the switch is noted as a reversed switch, but not added to the final filtered set
                            Data(c).VresetStatus(j) = "Net 0 flip";
                        elseif abs(Perc) > Thresh && Perc < 0
                            k = k+1;
                            Data(c).VresetInd(k) = Data(c).VresetRawIndex(j);
                            Data(c).Vreset(k) = Data(c).VresetRaw(j);
                            Data(c).VresetStatus(j) = "Net SET flip";
                        elseif abs(Perc) > Thresh && Perc > 0
                            k = k+1;
                            Data(c).VresetInd(k) = Data(c).VresetRawIndex(j);
                            Data(c).Vreset(k) = Data(c).VresetRaw(j);
                            Data(c).VresetStatus(j) = "Net RESET flip";
                        else
                            Data(c).VresetStatus(j) = "Error";
                        end
                    elseif Data(c).VresetDiff(j) < 0
                        Data(c).VresetStatus(j) = "Negative";
                    else
                        Data(c).VresetStatus(j) = "Error";
                    end
                end
                
                Data(c).R100R = min(interp1(V2, R2, 0.1, 'nearest'), interp1(V2, R2, -0.1, 'nearest'));
                Data(c).R100ratio = Data(c).R100F./Data(c).R100R;
                if Data(c).R100ratio > RatioMax %setting new Max ratio (and all associated parameters) if this particular scan shows a higher ratio, also replacing previous max with the new one
                    RatioMax2 = RatioMax;
                    Rmax2 = Rmax;
                    Rmin2 = Rmin;
                    RatioMax = Data(c).R100ratio;
                    Rmax = Data(c).R100F;
                    Rmin = Data(c).R100R;
                elseif (Data(c).R100ratio < RatioMax) && (Data(c).R100ratio > RatioMax2) % for finding 2nd highest on/off ratio.
                    RatioMax2 = Data(c).R100ratio;
                    Rmax2 = Data(c).R100F;
                    Rmin2 = Data(c).R100R;
                else
                end
            else
                Data(c).SweepType = "Single";
            end
        Data(c).Vmin = min(Data(c).Vout); %records minimum voltage
        Data(c).Vmax = max(Data(c).Vout);   %records maximum voltage
        end
    case 'No'
    otherwise
end

%The following loops take the first listed and detected SET/RESET voltages and put them in an array
k = 0;
l = 0;
VsetAgg = [];
VresetAgg = [];
for c = 1:size(Data, 2)
    if isfield(Data, 'Vset') == true
        if isempty(Data(c).Vset) == false
            k = k+1;
            VsetAgg(1, k) = Data(c).Vset(1);
        else
        end
    else
        VsetAgg = 0;
    end
    if isfield(Data, 'Vreset') == true
        if isempty(Data(c).Vreset) == false
            l = l+1;
            VresetAgg(1, l) = Data(c).Vreset(1);
        else
        end
    else
        VresetAgg = 0;
    end
end

for c = 1:size(Data, 2)
    if Data(c).R100F > R0
        R0 = Data(c).R100F;
    else
    end
end
        
%Now the software analyzes the array finding the average SET/RESET voltages and their standard errors
VsetAvg = mean(VsetAgg);
VsetErr = std(VsetAgg)/sqrt(size(VsetAgg, 1));
VresetAvg = mean(VresetAgg);
VresetErr = std(VresetAgg)/sqrt(size(VresetAgg, 1));

Values = [RatioMax Rmax Rmin RatioMax2 Rmax2 Rmin2 VsetAvg VsetErr VresetAvg VresetErr R0];%makes easily copy pasteable results of max ratio, resistance (high and low at the max ratio point), and starting resistance

%Past here is just figure plotting that can be fiddled with to see if the software potentially misread SET/RESET voltages or to see weird errors in the data.
if FigurePlot == 1
    %memristorplot(Data, save)
%{
    %    for x = 1:numel(Data)
%         savefileroot = convertStringsToChars(Data(x).Filepath);
%         savefileroot = savefileroot(1:end-4);
%     
%         tosave7 = figure;
%         ax = gca;
%         hold on
%         plot(Data(x).Vout, Data(x).I, '-k');
%         xlabel('Voltage (V)');
%         ylabel('Current (A)');
%        % axis([-4.1 2.6 -8E-3 8E-3]);
%         yl = ylim(ax);
%         axis(ax, 'tight')
%         ylim(ax, yl)
%         xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.05)
%         box on
%         set(gcf, 'Position', [100, 100, 600, 450])
%         set(gca, 'FontSize', 18)
%         hold off
%         if save == 1
%             saveas(tosave7,strcat(savefileroot,'_current.png'));
%             saveas(tosave7,strcat(savefileroot,'_current.fig'));
%         end
%         close
%         
%         tosave8 = figure;
%         ax = gca;
%         hold on
%         plot(Data(x).Vout, Data(x).Iabs, '-k');
%         xlabel('Voltage (V)');
%         ylabel('Current (A)');
% %         axis([-4.1 2.6 1.1E-8 1.4E-2]);
%         yl = ylim(ax);
%         axis(ax, 'tight')
%         ylim(ax, yl)
%         xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.05)
%         box on
%         set(gca, 'YScale', 'log');
%         set(gcf, 'Position', [200, 100, 600, 450])
%         set(gca, 'FontSize', 18)
%         hold off
%         if save == 1
%             saveas(tosave8,strcat(savefileroot,'_LogCurrent.png'));
%             saveas(tosave8,strcat(savefileroot,'_LogCurrent.fig'));
%         end
%         close
% 
%         tosave9 = figure;
%         ax = gca;
%         hold on
%         plot(Data(x).Vout, Data(x).R, '-k');
%         xlabel('Voltage (V)');
%         ylabel('Resistance (\Omega)');
%         %axis([-4 2.5 2E2 inf]);
%         yl = ylim(ax);
%         axis(ax, 'tight')
%         ylim(ax, yl)
%         xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.05)
%         box on
%         set(gca, 'YScale', 'log');
%         set(gcf, 'Position', [300, 100, 600, 450])
%         set(gca, 'FontSize', 18)
%         hold off
%         if save == 1
%             saveas(tosave9,strcat(savefileroot,'_Resistance.png'));
%             saveas(tosave9,strcat(savefileroot,'_Resistance.fig'));
%         end
%         close
%         
%          tosave9 = figure;
%         ax = gca;
%         hold on
%         plot(Data(x).Vout, Data(x).R, '-k');
%         xlabel('Voltage (V)');
%         ylabel('Resistance (\Omega)');
%         %axis([-4 2.5 2E2 inf]);
%         yl = ylim(ax);
%         axis(ax, 'tight')
%         ylim(ax, yl)
%         xlim(ax, xlim(ax) + [-1,1]*range(xlim(ax)).* 0.05)
%         box on
% %         set(gca, 'YScale', 'log');
%         set(gcf, 'Position', [300, 100, 600, 450])
%         set(gca, 'FontSize', 18)
%         hold off
%         if save == 1
%             saveas(tosave9,strcat(savefileroot,'_ResistanceLinear.png'));
%             saveas(tosave9,strcat(savefileroot,'_ResistanceLinear.fig'));
%         end
%         close
%      
%    end
%}   
   %writing the data taken on a single memristor device to a single file 
   excelfilemaker(path,Data)

   
   
%Plotting vset/reset and other info.
%{
    figure(10)
    hold on
    plot(Data(3).Vout, Data(3).R, '-k', Data(5).Vout, Data(5).R, '-r');
    xlabel('V (V)');
    ylabel('Resistance (Ohms)');
    axis([-4 2.5 2E2 inf]);
    box on
    legend('First Switch', 'Second Switch', 'Location', 'northeast');
    set(gca, 'YScale', 'log');
    set(gcf, 'Position', [300, 100, 600, 450])
    set(gca, 'FontSize', 18)
    hold off
    
    figure(11)
    hold on
    plot(Data(3).Vout, Data(3).Iabs, '-k', Data(5).Vout, Data(5).Iabs, '-r');
    xlabel('V (V)');
    ylabel('Current (A)');
    axis([-4.1 2.6 1.1E-8 1.4E-2]);
    box on
    set(gca, 'YScale', 'log');
    legend('First Switch', 'Second Switch', 'Location', 'northeast');
    set(gca, 'YScale', 'log');
    set(gcf, 'Position', [300, 100, 600, 450])
    set(gca, 'FontSize', 18)
    hold off
    
   
    figure(10)
    hold on
    for c = i1:i2
        plotC = plot(Data(c).Vout, Data(c).Iabs, '-', 'color',[0.6, 0.6, 0.6])
    end
    plotA = plot(Data(16).Vout, Data(16).Iabs, '-b', 'LineWidth', 2)
    legend([plotA plotC], 'Last On/Off', 'Other Sweeps', 'Location', 'southwest')
    xlabel('V (V)');
    ylabel('Current(A)');
    axis([-4.2 2.8 1E-8 9E-3]);
    yticks([1E-8 1E-7 1E-6 1E-5 1E-4 1E-3])
    box on
    set(gca, 'YScale', 'log');
    set(gcf, 'Position', [400, 100, 600, 450])
    set(gca, 'FontSize', 18)
    hold off
    
    
    for c = i1:i2
        HRS(c-i1+1) = Data(c).R100F;
        LRS(c-i1+1) = Data(c).R100R;
    end
    figure(11)
    hold on
    plot(HRS, '.b', 'MarkerSize', 14)
    plot(LRS, '.r', 'MarkerSize', 14)
    yline(min(HRS), '--b')
    yline(max(LRS), '--r')
    legend('HRS', 'LRS', 'Location', 'southwest')
    xlabel('Switching #');
    ylabel('Resistance (\Omega)');
    axis([-1 65 2E2 9E5]);
    box on
    set(gca, 'YScale', 'log');
    set(gcf, 'Position', [500, 100, 600, 450])
    set(gca, 'FontSize', 18)
    hold off
    
    for c = i1:i2 %selecting which set and reset voltages to make a vector out of, this was done manually
        VSet(c-2) = Data(c).VSet;
        VReset(c-2)=Data(c).VReset;
    end
    [PSet, VPSet] = ecdf(VSet);
    [PReset, VPReset] = ecdf(VReset);
    for c = 1:size(PReset)  %For inverting the cumulative probability distribution
        PReset(c) = 1-PReset(c);
    end
    figure(12)
    hold on
    plot(VPSet, PSet, '.r', VPReset, PReset, '.b', 'MarkerSize', 14)
    xlabel('Voltage (V)');
    ylabel('Cumulative Probability (%)')
    axis([-4 2.5 -0.01 1.01]);
    box on
    legend('V_{Set}', 'V_{Reset}');
    yticks([0 0.2 0.4 0.6 0.8 1]);
    yticklabels({'0', '20', '40', '60', '80', '100'});
    set(gcf, 'Position', [500, 100, 600, 450])
    set(gca, 'FontSize', 18)
    hold off
    %}
else
    
end
% disp('Press a key to end program')
% pause; %holds program here until a button is pressed
% close all; %closes all plots so you don't have to individually cross them out.