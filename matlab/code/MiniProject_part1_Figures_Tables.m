clc
clear all
close all

% load the datasets
disp('LOAD Data')

CurrentDir=pwd;
PathLoadData=[CurrentDir filesep 'Data'];
PathSaveFigures=[CurrentDir filesep 'Figures'];
PathSaveResults=[CurrentDir filesep 'Results'];
PathSaveTables=[CurrentDir filesep 'Tables'];

datasetName='Result_1.mat';

load(fullfile(PathSaveResults, datasetName), 'result');

disp('Data LOADED')
pause(0.5)

%% Parameters

Time_Window=2; % time window used to compute meam Vm, Vm SD and Vm FFT 

%% Make a table for average across cell class

CellClassMeans=[];
Mean_Table_Part1=[];

Cell_Class_List={'EXC', 'PV', 'VIP', 'SST'};

for i=1:4
    
    CellClassMeans.Cell_Class{i,1}=(cell2mat(Cell_Class_List(i))); 
    CellClassMeans.Firing_Rate_Mean(i,1)=mean(result.(cell2mat(Cell_Class_List(i))).Firing_Rate(:,1), 'omitnan')
    CellClassMeans.Firing_Rate_SD(i,1)=std(result.(cell2mat(Cell_Class_List(i))).Firing_Rate(:,1), 'omitnan');
    CellClassMeans.AP_Threshold_Mean(i,1)=mean(result.(cell2mat(Cell_Class_List(i))).AP_Threshold(:,1), 'omitnan');
    CellClassMeans.AP_Threshold_SD(i,1)=std(result.(cell2mat(Cell_Class_List(i))).AP_Threshold(:,1), 'omitnan');
    CellClassMeans.AP_Duration_Mean(i,1)=mean(result.(cell2mat(Cell_Class_List(i))).AP_Duration(:,1), 'omitnan');
    CellClassMeans.AP_Duration_SD(i,1)=std(result.(cell2mat(Cell_Class_List(i))).AP_Duration(:,1), 'omitnan');
    CellClassMeans.mean_Vm_Mean(i,1)=mean(result.(cell2mat(Cell_Class_List(i))).mean_Vm(:,1), 'omitnan');
    CellClassMeans.mean_Vm_SD(i,1)=std(result.(cell2mat(Cell_Class_List(i))).mean_Vm(:,1), 'omitnan');
    CellClassMeans.SD_Vm_Mean(i,1)=mean(result.(cell2mat(Cell_Class_List(i))).SD_Vm(:,1), 'omitnan');
    CellClassMeans.SD_Vm_SD(i,1)=std(result.(cell2mat(Cell_Class_List(i))).SD_Vm(:,1), 'omitnan');
    CellClassMeans.FFT_LF_Mean(i,1)=mean(result.(cell2mat(Cell_Class_List(i))).FFT_LF(:,1), 'omitnan');
    CellClassMeans.FFT_LF_SD(i,1)=std(result.(cell2mat(Cell_Class_List(i))).FFT_LF(:,1), 'omitnan');
    CellClassMeans.FFT_HF_Mean(i,1)=mean(result.(cell2mat(Cell_Class_List(i))).FFT_HF(:,1), 'omitnan');
    CellClassMeans.FFT_HF_SD(i,1)=std(result.(cell2mat(Cell_Class_List(i))).FFT_HF(:,1), 'omitnan');
end

Mean_Table_Part1=table(CellClassMeans.Cell_Class, CellClassMeans.Firing_Rate_Mean, CellClassMeans.Firing_Rate_SD, CellClassMeans.AP_Threshold_Mean, CellClassMeans.AP_Threshold_SD, ...
CellClassMeans.AP_Duration_Mean,CellClassMeans.AP_Duration_SD, CellClassMeans.mean_Vm_Mean,CellClassMeans.mean_Vm_SD,...
CellClassMeans.SD_Vm_Mean,CellClassMeans.SD_Vm_SD, CellClassMeans.FFT_LF_Mean,CellClassMeans.FFT_LF_SD, ...
CellClassMeans.FFT_HF_Mean, CellClassMeans.FFT_HF_SD,...
'VariableNames',{'Cell Class', 'Mean Firing Rate','± SD1', 'Mean AP Threshold','± SD2', 'Mean AP Duration','± SD3', 'Mean Vm','± SD4', ...
    'Mean Vm SD','± SD5', 'Mean LF FFT','± SD6', 'Mean HF FFT','± SD7',})
Expression=[PathSaveTables filesep 'Mean_Table_Part1.xls'];
writetable(Mean_Table_Part1, Expression)

%% Conpare mean Firing Rates across cells

clc

Tobeplotted_EXC=[];
Tobeplotted_PV=[];
Tobeplotted_VIP=[];
Tobeplotted_SST=[];

Tobeplotted_EXC(1:size(result.EXC.Firing_Rate,1),1)=1;
Tobeplotted_EXC(1:size(result.EXC.Firing_Rate,1),2)=result.EXC.Firing_Rate(:,1);

Tobeplotted_PV(1:size(result.PV.Firing_Rate,1),1)=2;
Tobeplotted_PV(1:size(result.PV.Firing_Rate,1),2)=result.PV.Firing_Rate(:,1);

Tobeplotted_VIP(1:size(result.VIP.Firing_Rate,1),1)=3;
Tobeplotted_VIP(1:size(result.VIP.Firing_Rate,1),2)=result.VIP.Firing_Rate(:,1);

Tobeplotted_SST(1:size(result.SST.Firing_Rate,1),1)=4;
Tobeplotted_SST(1:size(result.SST.Firing_Rate,1),2)=result.SST.Firing_Rate(:,1);


%%
figure

plot(Tobeplotted_EXC(:,1), Tobeplotted_EXC(:,2), 'O', 'Color', '[0.5 0.5 0.5]')
hold on
plot(Tobeplotted_PV(:,1), Tobeplotted_PV(:,2), 'O', 'Color', '[1 0.5 0.5]')
hold on
plot(Tobeplotted_VIP(:,1), Tobeplotted_VIP(:,2), 'O', 'Color', '[0.5 0.5 1]')
hold on
plot(Tobeplotted_SST(:,1), Tobeplotted_SST(:,2), 'O', 'Color', '[1 0.8 0.3]')
hold on
errorbar([1.2 2.2 3.2 4.2], CellClassMeans.Firing_Rate_Mean(:,1), CellClassMeans.Firing_Rate_SD(:,1),'O' , 'MarkerSize', 10,'Color', '[0 0 0]')

ax = gca;
ax.TickDir = 'out';
ax.XTick=[1, 2, 3, 4];
ax.XTickLabels={'EXC', 'PV', 'VIP', 'SST'};
xlim([0.5 5])
ylim([-1 60])
Graph_Title=['Mean Firing Rate'];
title(Graph_Title) % write the tittle of the graph
xlabel('Cell Class') % label the x axis
ylabel('FR (Hz)') % label the y axis

%% SAVE THE RESULT FIGURES

disp('Saving Figure')
pause(0.5)

Expression=[PathSaveFigures filesep '1_Mean_Firing_Rate'];

print('-painters', '-depsc', Expression)
print('-painters', '-djpeg', Expression)

disp('DONE')
pause(0.5)

%% Conpare mean AP duration across cells

clc

Tobeplotted_EXC=[];
Tobeplotted_PV=[];
Tobeplotted_VIP=[];
Tobeplotted_SST=[];

Tobeplotted_EXC(1:size(result.EXC.AP_Duration,1),1)=1;
Tobeplotted_EXC(1:size(result.EXC.AP_Duration,1),2)=result.EXC.AP_Duration(:,1);

Tobeplotted_PV(1:size(result.PV.AP_Duration,1),1)=2;
Tobeplotted_PV(1:size(result.PV.AP_Duration,1),2)=result.PV.AP_Duration(:,1);

Tobeplotted_VIP(1:size(result.VIP.AP_Duration,1),1)=3;
Tobeplotted_VIP(1:size(result.VIP.AP_Duration,1),2)=result.VIP.AP_Duration(:,1);

Tobeplotted_SST(1:size(result.SST.AP_Duration,1),1)=4;
Tobeplotted_SST(1:size(result.SST.AP_Duration,1),2)=result.SST.AP_Duration(:,1);


%%
figure

plot(Tobeplotted_EXC(:,1), Tobeplotted_EXC(:,2), 'O', 'Color', '[0.5 0.5 0.5]')
hold on
plot(Tobeplotted_PV(:,1), Tobeplotted_PV(:,2), 'O', 'Color', '[1 0.5 0.5]')
hold on
plot(Tobeplotted_VIP(:,1), Tobeplotted_VIP(:,2), 'O', 'Color', '[0.5 0.5 1]')
hold on
plot(Tobeplotted_SST(:,1), Tobeplotted_SST(:,2), 'O', 'Color', '[1 0.8 0.3]')
hold on
errorbar([1.2 2.2 3.2 4.2], CellClassMeans.AP_Duration_Mean(:,1), CellClassMeans.AP_Duration_SD(:,1),'O' , 'MarkerSize', 10,'Color', '[0 0 0]')

ax = gca;
ax.TickDir = 'out';
ax.XTick=[1, 2, 3, 4];
ax.XTickLabels={'EXC', 'PV', 'VIP', 'SST'};
xlim([0.5 5])
ylim([0 3])
Graph_Title=['Mean AP duration'];
title(Graph_Title) % write the tittle of the graph
xlabel('Cell Class') % label the x axis
ylabel('Duration (ms)') % label the y axis

%% SAVE THE RESULT FIGURES

disp('Saving Figure')
pause(0.5)

Expression=[PathSaveFigures filesep '2_Mean_AP_Duration'];

print('-painters', '-depsc', Expression)
print('-painters', '-djpeg', Expression)

disp('DONE')
pause(0.5)

%% Conpare mean Firing rate vs AP duration across cells

clc

Tobeplotted_EXC=[];
Tobeplotted_PV=[];
Tobeplotted_VIP=[];
Tobeplotted_SST=[];

Tobeplotted_EXC(1:size(result.EXC.Firing_Rate,1),1)=1;
Tobeplotted_EXC(1:size(result.EXC.Firing_Rate,1),2)=log10(result.EXC.Firing_Rate(:,1));
Tobeplotted_EXC(1:size(result.EXC.Firing_Rate,1),3)=result.EXC.AP_Duration(:,1);

Tobeplotted_PV(1:size(result.PV.Firing_Rate,1),1)=2;
Tobeplotted_PV(1:size(result.PV.Firing_Rate,1),2)=log10(result.PV.Firing_Rate(:,1));
Tobeplotted_PV(1:size(result.PV.Firing_Rate,1),3)=result.PV.AP_Duration(:,1);

Tobeplotted_VIP(1:size(result.VIP.Firing_Rate,1),1)=3;
Tobeplotted_VIP(1:size(result.VIP.Firing_Rate,1),2)=log10(result.VIP.Firing_Rate(:,1));
Tobeplotted_VIP(1:size(result.VIP.Firing_Rate,1),3)=result.VIP.AP_Duration(:,1);

Tobeplotted_SST(1:size(result.SST.Firing_Rate,1),1)=4;
Tobeplotted_SST(1:size(result.SST.Firing_Rate,1),2)=log10(result.SST.Firing_Rate(:,1));
Tobeplotted_SST(1:size(result.SST.Firing_Rate,1),3)=result.SST.AP_Duration(:,1);

%% Compute Pearson correlation between Firing rate, AP threshold, mean Vm and Vm SD
 
Temp=[];
All_Cells=[];

Temp=vertcat(Tobeplotted_EXC, Tobeplotted_PV, Tobeplotted_VIP, Tobeplotted_SST);
Temp(Temp==-Inf)=NaN;
All_Cells=rmmissing(Temp);

[Corr_R Corr_P]=corrcoef(All_Cells)

%%
fx=[];
fx=fit(All_Cells(:,3), All_Cells(:,2), 'poly1');


figure

plot(Tobeplotted_EXC(:,3), Tobeplotted_EXC(:,2), 'O', 'Color', '[0 0 0]')
hold on
plot(Tobeplotted_PV(:,3), Tobeplotted_PV(:,2), 'O', 'Color', '[1 0 0]')
hold on
plot(Tobeplotted_VIP(:,3), Tobeplotted_VIP(:,2), 'O', 'Color', '[0 0 1]')
hold on
plot(Tobeplotted_SST(:,3), Tobeplotted_SST(:,2), 'O', 'Color', '[1 0.5 0]')
hold on
plot(fx)

ax = gca;
ax.TickDir = 'out';
xlim([0 2.5])
ylim([-2.5 2.5])
Graph_Title=['Mean Firing rate vs AP duration'];
title(Graph_Title) % write the tittle of the graph
xlabel('AP Duration (ms)') % label the x axis
ylabel('Firing Rate (Hz)') % label the y axis
r=Corr_R(2, 3);
p = Corr_P(2, 3);
Expression=['r= ', num2str(r)];
text(ax,1.5,-1.8,Expression, 'FontSize',8)
Expression=['p= ', num2str(p)];
text(ax,1.5,-2,Expression, 'FontSize',8)

%% SAVE THE RESULT FIGURES

disp('Saving Figure')
pause(0.5)

Expression=[PathSaveFigures filesep '3_Mean_FRvsAP_Duration'];

print('-painters', '-depsc', Expression)
print('-painters', '-djpeg', Expression)

disp('DONE')
pause(0.5)

%% Conpare mean Vm across cells

clc

Tobeplotted_EXC=[];
Tobeplotted_PV=[];
Tobeplotted_VIP=[];
Tobeplotted_SST=[];

Tobeplotted_EXC(1:size(result.EXC.mean_Vm,1),1)=1;
Tobeplotted_EXC(1:size(result.EXC.mean_Vm,1),2)=result.EXC.mean_Vm(:,1);

Tobeplotted_PV(1:size(result.PV.mean_Vm,1),1)=2;
Tobeplotted_PV(1:size(result.PV.mean_Vm,1),2)=result.PV.mean_Vm(:,1);

Tobeplotted_VIP(1:size(result.VIP.mean_Vm,1),1)=3;
Tobeplotted_VIP(1:size(result.VIP.mean_Vm,1),2)=result.VIP.mean_Vm(:,1);

Tobeplotted_SST(1:size(result.SST.mean_Vm,1),1)=4;
Tobeplotted_SST(1:size(result.SST.mean_Vm,1),2)=result.SST.mean_Vm(:,1);


%%
figure

plot(Tobeplotted_EXC(:,1), Tobeplotted_EXC(:,2), 'O', 'Color', '[0.5 0.5 0.5]')
hold on
plot(Tobeplotted_PV(:,1), Tobeplotted_PV(:,2), 'O', 'Color', '[1 0.5 0.5]')
hold on
plot(Tobeplotted_VIP(:,1), Tobeplotted_VIP(:,2), 'O', 'Color', '[0.5 0.5 1]')
hold on
plot(Tobeplotted_SST(:,1), Tobeplotted_SST(:,2), 'O', 'Color', '[1 0.8 0.3]')
hold on
errorbar([1.2 2.2 3.2 4.2], CellClassMeans.mean_Vm_Mean(:,1), CellClassMeans.mean_Vm_SD(:,1),'O' , 'MarkerSize', 10,'Color', '[0 0 0]')

ax = gca;
ax.TickDir = 'out';
ax.XTick=[1, 2, 3, 4];
ax.XTickLabels={'EXC', 'PV', 'VIP', 'SST'};
xlim([0.5 5])
ylim([-0.075 -0.035])
Graph_Title=['Mean Vm'];
title(Graph_Title) % write the tittle of the graph
xlabel('Cell Class') % label the x axis
ylabel('Vm (mV)') % label the y axis

%% SAVE THE RESULT FIGURES

disp('Saving Figure')
pause(0.5)

Expression=[PathSaveFigures filesep '4_Mean_Vm'];

print('-painters', '-depsc', Expression)
print('-painters', '-djpeg', Expression)

disp('DONE')
pause(0.5)

%% Conpare mean Vm SD across cells

clc

Tobeplotted_EXC=[];
Tobeplotted_PV=[];
Tobeplotted_VIP=[];
Tobeplotted_SST=[];

Tobeplotted_EXC(1:size(result.EXC.SD_Vm,1),1)=1;
Tobeplotted_EXC(1:size(result.EXC.SD_Vm,1),2)=result.EXC.SD_Vm(:,1);

Tobeplotted_PV(1:size(result.PV.SD_Vm,1),1)=2;
Tobeplotted_PV(1:size(result.PV.SD_Vm,1),2)=result.PV.SD_Vm(:,1);

Tobeplotted_VIP(1:size(result.VIP.SD_Vm,1),1)=3;
Tobeplotted_VIP(1:size(result.VIP.SD_Vm,1),2)=result.VIP.SD_Vm(:,1);

Tobeplotted_SST(1:size(result.SST.SD_Vm,1),1)=4;
Tobeplotted_SST(1:size(result.SST.SD_Vm,1),2)=result.SST.SD_Vm(:,1);


%%
figure

plot(Tobeplotted_EXC(:,1), Tobeplotted_EXC(:,2), 'O', 'Color', '[0.5 0.5 0.5]')
hold on
plot(Tobeplotted_PV(:,1), Tobeplotted_PV(:,2), 'O', 'Color', '[1 0.5 0.5]')
hold on
plot(Tobeplotted_VIP(:,1), Tobeplotted_VIP(:,2), 'O', 'Color', '[0.5 0.5 1]')
hold on
plot(Tobeplotted_SST(:,1), Tobeplotted_SST(:,2), 'O', 'Color', '[1 0.8 0.3]')
hold on
errorbar([1.2 2.2 3.2 4.2], CellClassMeans.SD_Vm_Mean(:,1), CellClassMeans.SD_Vm_SD(:,1),'O' , 'MarkerSize', 10,'Color', '[0 0 0]')

ax = gca;
ax.TickDir = 'out';
ax.XTick=[1, 2, 3, 4];
ax.XTickLabels={'EXC', 'PV', 'VIP', 'SST'};
xlim([0.5 5])
ylim([0 0.012])
Graph_Title=['Mean Vm SD'];
title(Graph_Title) % write the tittle of the graph
xlabel('Cell Class') % label the x axis
ylabel('Vm SD (V)') % label the y axis

%% SAVE THE RESULT FIGURES

disp('Saving Figure')
pause(0.5)

Expression=[PathSaveFigures filesep '5_Mean_Vm_SD'];

print('-painters', '-depsc', Expression)
print('-painters', '-djpeg', Expression)

disp('DONE')
pause(0.5)

%% Plot GRD Average FFT for each cell class

SR_Vm=result.EXC.SR_Vm(1,1);
color=[0 0 0; 1 0 0; 0 0 1; 1 0.5 0];

figure

for i=1:4
    FFT_Mtrx=[];
    FFT_Mean=[];
    FFT_sem=[];
    for c=1:size(result.(cell2mat(Cell_Class_List(i))).FFT,1)
    
        FFT_Mtrx(c,:)=result.(cell2mat(Cell_Class_List(i))).FFT{c,1}';
   
    end
        
    Function_Plot_Mean_FFT(FFT_Mtrx, SR_Vm, color(i,:))
    
    eval(['GRD_AVG_FFT_' cell2mat(Cell_Class_List(i)) '_Mean=FFT_Mean;']);
    eval(['GRD_AVG_FFT_' cell2mat(Cell_Class_List(i)) '_sem=FFT_sem;']);
    
end

ax = gca;
ax.TickDir = 'out';
xlim([0.5 100])
ylim([0 0.0015])
Graph_Title=['GRD AVG FFTs'];
title(Graph_Title) % write the tittle of the graph
xlabel('Frequency (Hz)') % label the x axis
ylabel('Amplitude (V)') % label the y axis
%% SAVE THE RESULT FIGURES

disp('Saving Figure')
pause(0.5)

Expression=[PathSaveFigures filesep '6_GRD_AVG_FFTs'];

print('-painters', '-depsc', Expression)
print('-painters', '-djpeg', Expression)

disp('DONE')
pause(0.5)



%% Conpare mean LF FFT across cells

clc

Tobeplotted_EXC=[];
Tobeplotted_PV=[];
Tobeplotted_VIP=[];
Tobeplotted_SST=[];

Tobeplotted_EXC(1:size(result.EXC.FFT_LF,1),1)=1;
Tobeplotted_EXC(1:size(result.EXC.FFT_LF,1),2)=result.EXC.FFT_LF(:,1);

Tobeplotted_PV(1:size(result.PV.FFT_LF,1),1)=2;
Tobeplotted_PV(1:size(result.PV.FFT_LF,1),2)=result.PV.FFT_LF(:,1);

Tobeplotted_VIP(1:size(result.VIP.FFT_LF,1),1)=3;
Tobeplotted_VIP(1:size(result.VIP.FFT_LF,1),2)=result.VIP.FFT_LF(:,1);

Tobeplotted_SST(1:size(result.SST.FFT_LF,1),1)=4;
Tobeplotted_SST(1:size(result.SST.FFT_LF,1),2)=result.SST.FFT_LF(:,1);


%%
figure

plot(Tobeplotted_EXC(:,1), Tobeplotted_EXC(:,2), 'O', 'Color', '[0.5 0.5 0.5]')
hold on
plot(Tobeplotted_PV(:,1), Tobeplotted_PV(:,2), 'O', 'Color', '[1 0.5 0.5]')
hold on
plot(Tobeplotted_VIP(:,1), Tobeplotted_VIP(:,2), 'O', 'Color', '[0.5 0.5 1]')
hold on
plot(Tobeplotted_SST(:,1), Tobeplotted_SST(:,2), 'O', 'Color', '[1 0.8 0.3]')
hold on
errorbar([1.2 2.2 3.2 4.2], CellClassMeans.FFT_LF_Mean(:,1), CellClassMeans.FFT_LF_SD(:,1),'O' , 'MarkerSize', 10,'Color', '[0 0 0]')

ax = gca;
ax.TickDir = 'out';
ax.XTick=[1, 2, 3, 4];
ax.XTickLabels={'EXC', 'PV', 'VIP', 'SST'};
xlim([0.5 5])
ylim([0 0.0015])
Graph_Title=['Mean LF FFT'];
title(Graph_Title) % write the tittle of the graph
xlabel('Cell Class') % label the x axis
ylabel('Amplitude (V)') % label the y axis

%% SAVE THE RESULT FIGURES

disp('Saving Figure')
pause(0.5)

Expression=[PathSaveFigures filesep '7_Mean_FFT_LF'];

print('-painters', '-depsc', Expression)
print('-painters', '-djpeg', Expression)

disp('DONE')
pause(0.5)

%% Conpare mean Firing rate vs AP duration across cells

clc

Tobeplotted_EXC=[];
Tobeplotted_PV=[];
Tobeplotted_VIP=[];
Tobeplotted_SST=[];

Tobeplotted_EXC(1:size(result.EXC.Firing_Rate,1),1)=log10(result.EXC.Firing_Rate(:,1));
Tobeplotted_EXC(1:size(result.EXC.AP_Threshold,1),2)=result.EXC.AP_Threshold(:,1);
Tobeplotted_EXC(1:size(result.EXC.mean_Vm,1),3)=result.EXC.mean_Vm(:,1);
Tobeplotted_EXC(1:size(result.EXC.SD_Vm,1),4)=result.EXC.SD_Vm(:,1);

Tobeplotted_PV(1:size(result.PV.Firing_Rate,1),1)=log10(result.PV.Firing_Rate(:,1));
Tobeplotted_PV(1:size(result.PV.AP_Threshold,1),2)=result.PV.AP_Threshold(:,1);
Tobeplotted_PV(1:size(result.PV.mean_Vm,1),3)=result.PV.mean_Vm(:,1);
Tobeplotted_PV(1:size(result.PV.SD_Vm,1),4)=result.PV.SD_Vm(:,1);

Tobeplotted_VIP(1:size(result.VIP.Firing_Rate,1),1)=log10(result.VIP.Firing_Rate(:,1));
Tobeplotted_VIP(1:size(result.VIP.Firing_Rate,1),2)=result.VIP.AP_Threshold(:,1);
Tobeplotted_VIP(1:size(result.VIP.mean_Vm,1),3)=result.VIP.mean_Vm(:,1);
Tobeplotted_VIP(1:size(result.VIP.SD_Vm,1),4)=result.VIP.SD_Vm(:,1);

Tobeplotted_SST(1:size(result.SST.Firing_Rate,1),1)=log10(result.SST.Firing_Rate(:,1));
Tobeplotted_SST(1:size(result.SST.AP_Threshold,1),2)=result.SST.AP_Threshold(:,1);
Tobeplotted_SST(1:size(result.SST.mean_Vm,1),3)=result.SST.mean_Vm(:,1);
Tobeplotted_SST(1:size(result.SST.SD_Vm,1),4)=result.SST.SD_Vm(:,1);
%% Compute Pearson correlation between Firing rate, AP threshold, mean Vm and Vm SD
 
Temp=[];
All_Cells=[];

Temp=vertcat(Tobeplotted_EXC, Tobeplotted_PV, Tobeplotted_VIP, Tobeplotted_SST);
Temp(Temp==-Inf)=NaN;
All_Cells=rmmissing(Temp);

[Corr_R Corr_P]=corrcoef(All_Cells)

%%
figure('Position', [100 300 1200 300])

% plot Firing rate vs AP Threshold

fx=fit(All_Cells(:,2), All_Cells(:,1), 'poly1');

subplot(1,3,1)
plot(Tobeplotted_EXC(:,2), Tobeplotted_EXC(:,1), '^', 'Color', '[0 0 0]')
hold on
plot(Tobeplotted_PV(:,2), Tobeplotted_PV(:,1), 'O', 'Color', '[1 0 0]')
hold on
plot(Tobeplotted_VIP(:,2), Tobeplotted_VIP(:,1), 'O', 'Color', '[0 0 1]')
hold on
plot(Tobeplotted_SST(:,2), Tobeplotted_SST(:,1), 'O', 'Color', '[1 0.5 0]')
hold on
plot(fx)

ax = gca;
ax.TickDir = 'out';
xlim([-0.055 -0.025])
ylim([-3 3])
Graph_Title=['Mean Firing rate vs AP threshold'];
title(Graph_Title) % write the tittle of the graph
xlabel('AP Threshold (V)') % label the x axis
ylabel('log Firing Rate (Hz)') % label the y axis
r=Corr_R(1, 2);
p = Corr_P(1, 2);
Expression=['r= ', num2str(r)];
text(ax,-0.05,2.8,Expression, 'FontSize',8)
Expression=['p= ', num2str(p)];
text(ax,-0.05,2.3,Expression, 'FontSize',8)

% plot Firing rate vs Mean Vm

fx=fit(All_Cells(:,3), All_Cells(:,1), 'poly1');

subplot(1,3,2)
plot(Tobeplotted_EXC(:,3), Tobeplotted_EXC(:,1), '^', 'Color', '[0 0 0]')
hold on
plot(Tobeplotted_PV(:,3), Tobeplotted_PV(:,1), 'O', 'Color', '[1 0 0]')
hold on
plot(Tobeplotted_VIP(:,3), Tobeplotted_VIP(:,1), 'O', 'Color', '[0 0 1]')
hold on
plot(Tobeplotted_SST(:,3), Tobeplotted_SST(:,1), 'O', 'Color', '[1 0.5 0]')
hold on
plot(fx)

ax = gca;
ax.TickDir = 'out';
xlim([-0.075 -0.035])
ylim([-3 3])
Graph_Title=['Mean Firing rate vs Mean Vm'];
title(Graph_Title) % write the tittle of the graph
xlabel('Mean Vm (V)') % label the x axis
ylabel('log Firing Rate (Hz)') % label the y axis
r=Corr_R(1, 3);
p = Corr_P(1, 3);
Expression=['r= ', num2str(r)];
text(ax,-0.07,2.8,Expression, 'FontSize',8)
Expression=['p= ', num2str(p)];
text(ax,-0.07,2.3,Expression, 'FontSize',8)

% plot Firing rate vs Vm SD

fx=fit(All_Cells(:,4), All_Cells(:,1), 'poly1');
subplot(1,3,3)
plot(Tobeplotted_EXC(:,4), Tobeplotted_EXC(:,1), '^', 'Color', '[0 0 0]')
hold on
plot(Tobeplotted_PV(:,4), Tobeplotted_PV(:,1), 'O', 'Color', '[1 0 0]')
hold on
plot(Tobeplotted_VIP(:,4), Tobeplotted_VIP(:,1), 'O', 'Color', '[0 0 1]')
hold on
plot(Tobeplotted_SST(:,4), Tobeplotted_SST(:,1), 'O', 'Color', '[1 0.5 0]')
hold on
plot(fx)

ax = gca;
ax.TickDir = 'out';
xlim([0 0.011])
ylim([-3 3])
Graph_Title=['Mean Firing rate vs Vm SD'];
title(Graph_Title) % write the tittle of the graph
xlabel('Vm SD (V)') % label the x axis
ylabel('log Firing Rate (Hz)') % label the y axis
r=Corr_R(1, 4);
p = Corr_P(1, 4);
Expression=['r= ', num2str(r)];
text(ax,0.001,2.8,Expression, 'FontSize',8)
Expression=['p= ', num2str(p)];
text(ax,0.001,2.3,Expression, 'FontSize',8)


%% SAVE THE RESULT FIGURES

disp('Saving Figure')
pause(0.5)

Expression=[PathSaveFigures filesep '8_Mean_FR_Correlations'];

print('-painters', '-depsc', Expression)
print('-painters', '-djpeg', Expression)

disp('DONE')
pause(0.5)


