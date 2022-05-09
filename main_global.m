%% Global calibration of Hs and Tm from WAVERYS and ERA5 wave reanalyses

%% Import data
clc, clear, close all
topLevelFolder = 'FOLDER_CONTAINING_DATA';
cd(topLevelFolder);
addpath(strcat(topLevelFolder,'/Functions'));
load('Data_Faro_Buoy_WAVERYS_ERA5_1993_2019.mat');
Var_Na = ["Hs","Tm"];
Var_Un = ["m","s"];
Model = string({'WAVERYS';'ERA5'});
% Parameters for the global calibration
% WAVERYS) Hs_cal = a*Hs^b; Tm_cal = Tm*p1+p2;
% ERA5) Hs_cal = a*Hs^b; Tm_cal = Tm*p1+p2;
power.a.WAVERYS = 1.045;
power.b.WAVERYS = 1.021;
poly.p1.WAVERYS = 0.8997;
poly.p2.WAVERYS = 1.011;
power.a.ERA5 = 1.073;
power.b.ERA5 = 1.02;
poly.p1.ERA5 = 0.945;
poly.p2.ERA5 = 1.124;

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
%% Calibrating Hs and Tm with the global equations
for ii = 1:numel(Var_Na) %loop over the variables (1 Hs, 2 Tm)
    for i = 1:numel(Model) %loop over the models (1 WAVERYS, 2 ERA5)
        if ii == 1 % Global calibration of Hs
            %Power fit (y=ax^b)
            TT{ii,i}.Glofit = power.a.(Model(i)).*TT{ii,i}.Model.^power.b.(Model(i));
        elseif ii == 2 % Global calibration of Tm
            %Polinomial fit 1st order with intercept (y=ax+b)
            TT{ii,i}.Glofit = TT{ii,i}.Model.*poly.p1.(Model(i))+poly.p2.(Model(i));        
        end
    end
end
clear i ii

%%
%% Calculating skill metrics for original and calibrated reanalysis data versus buoy data
STAT = struct;
for ii = 1:numel(Var_Na) %loop over the variables (1 Hs, 2 Tm)
    for i = 1:numel(Model) %loop over the models (1 WAVERYS, 2 ERA5)
        STAT(i).Var(ii).cal = ValStat(TT{ii,i}.Buoy,TT{ii,i}.Glofit); 
        STAT(i).Var(ii).ori = ValStat(TT{ii,i}.Buoy,TT{ii,i}.Model); 
    end
end
clear ii i j
%%
%%
%% Plot of original and calibrated reanalysis data
for i = 1:numel(Model) %loop over the models (1 WAVERYS, 2 ERA5)
    figure;
    t = tiledlayout(1,4);
    t.TileSpacing = 'tight';
    set(gcf,'position',[0,0,1500,350]);
    title(t,strcat(' Buoy -  ',Model(i)));
    for ii = 1:numel(Var_Na) %loop over the variables (Hs 1, Tm 2)
        maxlim(ii)=max(max(max(TT{ii,1}.Buoy),max(TT{ii,1}.Model)),max(max(TT{ii,2}.Buoy),max(TT{ii,2}.Model)));
        % Create vectors with percentiles
        Prctiles = zeros(99,2);
        for k = 1:99
            Prctiles(k,1) = prctile(TT{ii,i}.Buoy,k); % Buoy
            Prctiles(k,2) = prctile(TT{ii,i}.Model,k); % Model
        end
        T_Prctiles = array2table(Prctiles,'VariableNames',{'Buoy','Model'});
        
        nexttile;
        h1=binscatter(TT{ii,i}.Buoy,TT{ii,i}.Model,[250 250]);
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
        title(sprintf('Uncalibrated Model - %s',string(Model(i))));
        line_leg = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        string_leg=sprintf('R=%0.2f\n RMSE=%0.2f\n Bias=%0.2f\n StDev=%0.2f\n SI=%0.2f',...
        STAT(i).Var(ii).ori.R,STAT(i).Var(ii).ori.RMSE,STAT(i).Var(ii).ori.Bias,STAT(i).Var(ii).ori.StDev,STAT(i).Var(ii).ori.SI);
        h_leg=legend(line_leg,string_leg,'Location','northwest');
        set(h_leg,'visible','off')
        annotation('textbox',(get(h_leg,'position')+[0 0 0 0]), 'String', string_leg,'FitBoxToText','on')
        clear fit line_leg string_leg h_leg h1 bs T_Prctiles Prctiles
        
        % Create vectors with percentiles
        Prctiles = zeros(99,2);
        for k = 1:99
            Prctiles(k,1) = prctile(TT{ii,i}.Buoy,k); % Buoy
            Prctiles(k,2) = prctile(TT{ii,i}.Glofit,k); % Model
        end
        T_Prctiles = array2table(Prctiles,'VariableNames',{'Buoy','Model'});
        nexttile;
        h1=binscatter(TT{ii,i}.Buoy,TT{ii,i}.Glofit,[250 250]);
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
        title(sprintf('Calibrated Model - %s',string(Model(i))));
        line_leg = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        string_leg=sprintf('R=%0.2f\n RMSE=%0.2f\n Bias=%0.2f\n StDev=%0.2f\n SI=%0.2f',...
        STAT(i).Var(ii).cal.R,STAT(i).Var(ii).cal.RMSE,STAT(i).Var(ii).cal.Bias,STAT(i).Var(ii).cal.StDev,STAT(i).Var(ii).cal.SI);
        h_leg=legend(line_leg,string_leg,'Location','northwest');
        set(h_leg,'visible','off')
        annotation('textbox',(get(h_leg,'position')+[0 0 0 0]), 'String', string_leg,'FitBoxToText','on')
        clear fit line_leg string_leg h_leg h1 bs T_Prctiles Prctiles
    end
    clear t
    %saveas(gcf, fullfile(PlotsFolder,strcat('Figure',num2str(i))));
end
 clear i ii j x1 h1