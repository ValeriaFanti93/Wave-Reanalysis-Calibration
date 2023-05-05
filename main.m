%% Calibration of WAVERYS and ERA5 wave reanalyses with buoy data
% The present script allows for validating and calibrating wave reanalysis
% data against buoy observations.

% The input data for the buoy, WAEVERYS and ERA5 are organized in a unique file with a timetable with:
% - time (yyyy-mm-dd HH:MM:SS)
% - significant wave height (Hs);
% - mean wave period (Tm);
% - peak wave period (Tp);
% - mean wave direction (Dir).

%% Import data
clc, clear, close all
% Starting folder
topLevelFolder = 'FOLDER_CONTAINING_DATA';
cd(topLevelFolder);
addpath(strcat(topLevelFolder,'/Functions'));
PlotsFolder = strcat(topLevelFolder,'/Plots');
file_name = 'Data_GL_TS_MO_42035_WAVERYS_ERA5_1993_2020.mat';
file_data = dir(file_name);
load(file_data.name);
Var_Na = ["Hs","Tm","Tp","Dir"];
Var_Un = ["m","s","s","°"];
Model = {'WAVERYS','ERA5'};
clear topLevelFolder file_name file_data 
%%
%% Load Timeseries
% Removing NaNs from buoy data and creating a timevector with same dates
% for buoy and reanalyses by retiming them.

for i = 1:numel(Var_Na)
    BUOY = rmmissing(timetable(DATA.BUOY.Time, DATA.BUOY.(Var_Na(i)))); 
    WAVERYS = retime(timetable(DATA.WAVERYS.Time,DATA.WAVERYS.(Var_Na(i))),BUOY.Time);
    ERA5 = retime(timetable(DATA.ERA5.Time,DATA.ERA5.(Var_Na(i))),BUOY.Time);
    TT{i,1} = timetable(BUOY.Time,BUOY.Var1,WAVERYS.Var1);
    % Remove NaNs because WAVERYS is 3 hourly
    TT{i,1} = rmmissing(TT{i,1}); 
    TT{i,1}.Properties.VariableNames = {'Buoy','Model'};
    TT{i,2} = timetable(BUOY.Time,BUOY.Var1,ERA5.Var1);
    TT{i,2}.Properties.VariableNames = {'Buoy','Model'};
end
 clear i BUOY WAVERYS ERA5
%%
%% Applying calibration to reanalysis data
Fit_name = {'NoFit','Linear','Polinomial1','Power1','Rotation'};
Fit_eq = {'y=y','y=ax','y=ax+b','y=ax^b','y=ROTx'};
Transf = struct;
for ii = 1:numel(Var_Na)-1 %loop over the variables (1 Hs, 2 Tm, 3 Tp)
     for i=1:numel(Model) %loop over the models (1 WAVERYS, 2 ERA5)
        TTModelCAL = CalReanalysis(TT{ii,i}.Buoy,TT{ii,i}.Model);
        for j = 1:width(TTModelCAL)
            fit = strcat('fit',num2str(j));
            TT{ii,i}.(fit) = TTModelCAL(:,j);
            clear fit
        end
        clear TTModelCAL
     end
end
clear i ii
%%
%% Calculating skill metrics for original and calibrated reanalysis data versus buoy data
STAT = struct;
for ii=1 : numel(Var_Na)-1 %loop over the variables (1 Hs, 2 Tm, 3 Tp)
    for i=1:numel(Model) %loop over the models (1 WAVERYS, 2 ERA5)
        for j=1:numel(Fit_name)
            fit = strcat('fit',num2str(j));
            STAT(i).Var(ii).fit(j) = ValStat(TT{ii,i}.Buoy,TT{ii,i}.(fit));
            clear fit
        end    
    end
end
clear ii i j gcf 
%%
%% Plot of original and calibrated reanalysis data
for i = 1:numel(Model) %loop over the models (1 WAVERYS, 2 ERA5)
    figure;
    t = tiledlayout(3,5);
    set(gcf,'position',[0,0,1500,950]);
    title(t,strcat(' Buoy -  ',Model(i)));
    for ii = 1:numel(Var_Na)-1 %loop over the variables (Hs 1, Tm 2, Tp 3)
        maxlim(ii)=max(max(max(TT{ii,1}.Buoy),max(TT{ii,1}.Model)),max(max(TT{ii,2}.Buoy),max(TT{ii,2}.Model)));
        for j = 1:length(Fit_name)
            fit = strcat('fit',num2str(j));
            % Create vectors with percentiles
            Prctiles = zeros(99,2);
            for k = 1:99
                Prctiles(k,1) = prctile(TT{ii,i}.Buoy,k); % Buoy
                Prctiles(k,2) = prctile(TT{ii,i}.(fit),k); % Model
            end
            T_Prctiles = array2table(Prctiles,'VariableNames',{'Buoy','Model'});
            
            nexttile;
            h1=binscatter(TT{ii,i}.Buoy,TT{ii,i}.(fit),[250 250]);
            set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            colormap(gca,'parula');
            bs=gca;
            bs.ColorScale = 'log';
            xlim([0 maxlim(ii)])
            ylim([0 maxlim(ii)])
            cline=refline(1,0);
            cline.Color='k';
            hold on 
            scatter(T_Prctiles.Buoy,T_Prctiles.Model,8,'filled','MarkerEdgeColor','r','MarkerFaceColor','r');
            hold off
            xlabel(strcat(Var_Na(ii),sprintf(' Buoy (%s)',Var_Un(ii))));
            ylabel(strcat(Var_Na(ii),sprintf('(%s)',Var_Un(ii))));
            if j == 1
                title(sprintf('Uncalibrated Model - %s',string(Fit_eq(j))));
            else
                title(sprintf('Calibrated Model - %s',string(Fit_eq(j))));
            end
            line_leg = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            string_leg=sprintf('R=%0.2f\n RMSE=%0.2f\n Bias=%0.2f\n StDev=%0.2f\n SI=%0.2f',...
            STAT(i).Var(ii).fit(j).R,STAT(i).Var(ii).fit(j).RMSE,STAT(i).Var(ii).fit(j).Bias,STAT(i).Var(ii).fit(j).StDev,STAT(i).Var(ii).fit(j).SI);
            h_leg=legend(line_leg,string_leg,'Location','northwest');
            set(h_leg,'visible','off')
            annotation('textbox',(get(h_leg,'position')+[0 0 0 0]), 'String', string_leg,'FitBoxToText','on')
            clear fit line_leg string_leg h_leg h1 bs T_Prctiles Prctiles
        end
    end
    clear t
    saveas(gcf, fullfile(PlotsFolder,strcat('Figure',num2str(i))));
end
 clear i ii j x1 h1
%%
%% Choosing the best fit for each variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To change. The calibrating equation here chosen are strictly related to the 
% example shown and must be changed when appying the code to a differnt buoy.
% ChooseFit equals 1 if no calibration is needed, 2 for a first degree polynomial
% calibration without intercept, 3 for a first degree polynomial calibration with
% intercept, 4 for a power function, 5 for a rotation around the mean.
ChooseFit(1,1) = 5;% Hs WAVERYS
ChooseFit(1,2) = 5;% Hs ERA5
ChooseFit(2,1) = 4;% Tm WAVERYS
ChooseFit(2,2) = 4;% Tm ERA5
ChooseFit(3,1) = 1;% Tp WAVERYS
ChooseFit(3,2) = 1;% Tp ERA5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saving calibrated variables in a structure
for i=1:numel(Model) %loop over the models (1 WAVERYS, 2 ERA5)
    for ii=1:numel(Var_Na)-1 %loop over the variables (1 Hs, 2 Tm, 3 Tp)
        fit = strcat('fit',num2str(ChooseFit(ii,i)));
        CAL{1,i}.(Var_Na(ii)) = timetable(TT{ii,i}.Time,TT{ii,i}.(fit));
        clear fit
    end
end
clear i ii
%%
%% Plot of best fit
for i=1:numel(Model) %loop over the models (1 WAVERYS, 2 ERA5)
    figure;
    t=tiledlayout(1,3);
    set(gcf,'position',[0,0,1200,400]);
    title(t,Model(i));
    for ii=1:numel(Var_Na)-1 %loop over the variables (Hs 1, Tm 2, Tp 3)
        nexttile;
        h1=binscatter(TT{ii,i}.Buoy,CAL{1,i}.(Var_Na(ii)).Var1,[250 250]);
        set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        colormap(gca,'parula');
        hold on
        bs=gca;
        bs.ColorScale = 'log';
        xlim([0 maxlim(ii)])
        ylim([0 maxlim(ii)])
        cline=refline(1,0);
        cline.Color='k';
        plot(prctile(TT{ii,i}.Buoy,1:99),prctile(CAL{1,i}.(Var_Na(ii)).Var1,1:99),'o','MarkerEdgeColor','r');
        hold off
        xlabel(strcat(Var_Na(ii),sprintf(' Buoy (%s)',Var_Un(ii))));
        ylabel(strcat(Var_Na(ii),strcat(' Model (',Var_Un(ii),')')));
        if ChooseFit(ii,i) == 1
            title(sprintf('Uncalibrated %s %s',string(Model(i)),string(Fit_eq(ChooseFit(ii,i))))); 
        else
            title(sprintf('Calibrated %s %s',string(Model(i)),string(Fit_eq(ChooseFit(ii,i)))));  
        end
        line_leg = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        string_leg=sprintf('R=%0.2f\n RMSE=%0.2f\n Bias=%0.2f\n StDev=%0.2f\n SI=%0.2f',...
        STAT(i).Var(ii).fit(ChooseFit(ii,i)).R,STAT(i).Var(ii).fit(ChooseFit(ii,i)).RMSE,STAT(i).Var(ii).fit(ChooseFit(ii,i)).Bias,...
        STAT(i).Var(ii).fit(ChooseFit(ii,i)).StDev,STAT(i).Var(ii).fit(ChooseFit(ii,i)).SI);
        h_leg=legend(line_leg,string_leg,'Location','southeast');
        set(h_leg,'visible','off')
        annotation('textbox',(get(h_leg,'position')+[0.01 0.03 0 0]), 'String', string_leg,'FitBoxToText','on')
        clear string_leg h_leg h1 bs
    end
   saveas(gcf, fullfile(PlotsFolder,strcat('Figure',num2str(i+2))));
   clear t
end
clear i ii
clear maxlim h1 bs line_leg h_leg string_leg
%%
%% Calibration of Directions
tiledlayout(1,2);
set(gcf,'position',[0,0,900,400]);
for i = 1:numel(Model)
    nexttile;
    h1 = histogram(TT{4,i}.Buoy);
    h1.FaceColor = [0 0 1];%[0 0.5 0.5];
    h1.BinWidth = 5;
    hold on
    h2 = histogram(TT{4,i}.Model);
    h2.FaceColor = [1 1 0];
    h2.BinWidth = 5;
    hold off
    xlabel('Wave direction (°)')
    ylabel('Occurrences')
    legend('Buoy',string(Model(i)),'Location','northwest');
end
saveas(gcf, fullfile(PlotsFolder,'Figure5'));
clear i h1 h2
%%
%% Choosing direction sectors for calibration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% To change. The directions sectors here chosen are strictly related to the 
% example shown and must be changed when appying the code to a differnt buoy.
dir1=0; dir2=180; dir3=360; dir4=0; dir5=0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:numel(Model) %loop over the models
    if dir1 == 0 && dir2 == 0 %No calibration needed
        [DirModelCal] = TT{4,i}.Model';
    elseif dir3 == 0
        [DirModelCal] = DirCal(TT{4,i}.Buoy,TT{4,i}.Model,[dir1,dir2]);   
    elseif dir4 == 0
        [DirModelCal] = DirCal(TT{4,i}.Buoy,TT{4,i}.Model,[dir1,dir2,dir3]);
    else
        [DirModelCal] = DirCal(TT{4,i}.Buoy,TT{4,i}.Model,[dir1,dir2,dir3,dir4]);
    end
% Saving calibrated directions in a timetable
CAL{1,i}.Dir = timetable(TT{4,i}.Time,DirModelCal');
clear DirModelCal
end
clear i dir1 dir2 dir3 dir4 dir5
%%
%% Plot of calibrated directions
tiledlayout(1,2);
set(gcf,'position',[0,0,900,400]);
for i = 1:numel(Model)
    nexttile;
    h1 = histogram(TT{4,i}.Buoy);
    h1.FaceColor = [0 0 1];%[0 0.5 0.5];
    h1.BinWidth = 5;
    hold on
    h2 = histogram(TT{4,i}.Model);
    h2.FaceColor = [1 1 0];
    h2.BinWidth = 5;
    h3 = histogram(CAL{1,i}.Dir.Var1);
    h3.FaceColor = [0 1 1];
    h3.BinWidth = 5;
    hold off
    xlabel('Wave direction (°)')
    ylabel('Occurrences')
    legend('Buoy',string(Model(i)),strcat('Calibrated',string(Model(i))),'Location','northwest');
end
saveas(gcf, fullfile(PlotsFolder,'Figure6'));
clear i h1 h2 h3
