clc
clear all
close all

if ~ exist('Results')
    
    mkdir('Results')
    addpath('Results')

end


if ~ exist('SavedFigures')
    
    mkdir('SavedFigures')
    addpath('SavedFigures')

end

if ~ exist('Tables')
    
    mkdir('Tables')
    addpath('Tables')

end


% load the datasets
disp('LOAD Data')

CurrentDir=pwd;
PathLoadData=[CurrentDir filesep 'Data'];
PathSaveFigures=[CurrentDir filesep 'SavedFigures'];
PathSaveResults=[CurrentDir filesep 'Results'];

datasetName='Data_Bio482.mat';

load(fullfile(PathLoadData, datasetName), 'Data');

disp('Data LOADED')
pause(0.5)


%% Parameters

Sweep_Type='free whisking'; % or 'active touch'
Cell_Types={'EXC', 'PV', 'VIP', 'SST'};
result=[];

TimeWindow=2; % Time window to analyze Vm (s)
FreqBandLim= [1 10 30 90]; % Low- and High-frequency Band limits (Hz)

%% Select for one Sweep type
dataSelect=[];
myFields = fieldnames(Data);
for thisField = 1:length(myFields)
    dataSelect.(cell2mat(myFields(thisField))) = Data.(cell2mat(myFields(thisField)))(strcmp(Data.Sweep_Type, Sweep_Type),:);
end

%% Print example Vm with AP peak and threshold

sweep=120;
disp(['Cell_ID = ', Data.Cell_ID{sweep}])

MembranePotential=[];
SR_Vm=[];
AP_Vm_Deriv_Thrs=[];
AP_Param=[];
Rec_Dur=[];

MembranePotential=Data.Sweep_MembranePotential{sweep,1};
SR_Vm=Data.Sweep_MembranePotential_SamplingRate(sweep,1);
AP_Vm_Deriv_Thrs=Data.Cell_APThreshold_Slope(1,1);
[AP_Param]=Function_Detect_APs(MembranePotential, SR_Vm, AP_Vm_Deriv_Thrs);
Rec_Dur=length(MembranePotential)/SR_Vm;
time_vect=linspace(0, Rec_Dur, length(MembranePotential));

figure('Position', [100 400 1400 400])
plot(time_vect, MembranePotential, 'Color', [0 0 0])
hold on
plot(AP_Param(:,1), AP_Param(:,2), 'o', 'Color', [1 0 0])
hold on
plot(AP_Param(:,3), AP_Param(:,4), 'o', 'Color', [0 0 1])
xlim([0 10])
xlabel('Time (s)')
ylabel('Vm (V)')
title(['cell ID = ' Data.Cell_ID{sweep}])
%% SAVE THE RESULT FIGURES

disp('Saving Figure')
pause(0.5)

Expression=[PathSaveFigures filesep '0_Example_Vm_AP_Thresholds'];

print('-painters', '-depsc', Expression)
print('-painters', '-djpeg', Expression)

disp('DONE')
pause(0.5)

%% Loop through Cell types

for tp=1:size(Cell_Types,2) % loop through the cell types
    
    % select one cell type
    
    Cell_Type=Cell_Types{tp};
    myFields = fieldnames(dataSelect);
    for thisField = 1:length(myFields)
        data1Type.(cell2mat(myFields(thisField))) = dataSelect.(cell2mat(myFields(thisField)))(strcmp(dataSelect.Cell_Type, Cell_Type),:);
    end
    
    Cell_List=unique(data1Type.Cell_ID);
    
    % Loop through Cells
    
    for c=1:size(Cell_List,1) % loop through the cells from 1 cell type
        
        % select 1 cell
        data1Cell=[];
        Cell_Name=Cell_List{c,1};
        myFields = fieldnames(data1Type);
        for thisField = 1:length(myFields)
            data1Cell.(cell2mat(myFields(thisField))) = data1Type.(cell2mat(myFields(thisField)))(strcmp(data1Type.Cell_ID, Cell_Name)==1, :);
        end
        
        % Fill result structure with metadata
        
        result.(cell2mat(Cell_Types(tp))).Cell_Name{c,1}=data1Cell.Cell_ID{1,1};
        result.(cell2mat(Cell_Types(tp))).Cell_Type{c,1}=data1Cell.Cell_Type{1,1};
        result.(cell2mat(Cell_Types(tp))).Cell_Depth(c,1)=data1Cell.Cell_Depth(1,1);
        result.(cell2mat(Cell_Types(tp))).Cell_Layer{c,1}=data1Cell.Cell_Layer{1,1};
        result.(cell2mat(Cell_Types(tp))).Cell_Anatomy{c,1}=data1Cell.Cell_Anatomy{1,1};
        result.(cell2mat(Cell_Types(tp))).Cell_Anatomy{c,2}=data1Cell.Cell_Anatomy{1,2};
        result.(cell2mat(Cell_Types(tp))).SR_Vm(c,1)=data1Cell.Sweep_MembranePotential_SamplingRate(1,1);
        
        % Initialyse variables
        Tot_Rec_Dur=0;
        Tot_Numb_APs=0;
        Tot_AP_Thrs=[];
        Tot_AP_Dur=[];
        Tot_Mean_Vm=[];
        Tot_SD_Vm=[];
        Tot_FFT=[];
        
        % Loop through sweeps for 1 cell
        
        for sweep=1:size(data1Cell.Sweep_Counter,1) % loop through the sweeps from 1 cell
            
            % get the data to be processed from the data structure
            
            MembranePotential=[];
            SR_Vm=[];
            AP_Vm_Deriv_Thrs=[];
            Rec_Dur=[];
            
            MembranePotential=data1Cell.Sweep_MembranePotential{sweep,1};
            SR_Vm=data1Cell.Sweep_MembranePotential_SamplingRate(sweep,1);
            AP_Vm_Deriv_Thrs=data1Cell.Cell_APThreshold_Slope(1,1);
            
            % update recording duration
            Rec_Dur=length(MembranePotential)/SR_Vm;
            Tot_Rec_Dur=Tot_Rec_Dur+Rec_Dur;
            
            % detect APs and compute AP parameters (1= AP Threshold times; 2= AP Threshols Vm; 3= AP peak times; 4= AP peak Vm; 5= AP amplitude)            
            
            [AP_Param]=Function_Detect_APs(MembranePotential, SR_Vm, AP_Vm_Deriv_Thrs);
            
            Tot_Numb_APs=Tot_Numb_APs+size(AP_Param,1);
            
            % Initialyse variables
                AP_Thrs_Times=[];
                AP_Thrs_Vm=[];
                AP_Peak_Times=[];
                AP_Peak_Vm=[];
                AP_Duration=[];
            
            % If APs have been detected => cut APs from Vm to obtain
            % subthreshold Vm
            
            if ~isempty(AP_Param)

                if min(AP_Param(:,6))<0

                    disp(['Error AP duration, Cell=', data1Cell.Cell_ID{1,1}])

                end

                Vm_Sub=[];
         
                AP_Thrs_Times=AP_Param(:,1);
                AP_Thrs_Vm=AP_Param(:,2);
                AP_Peak_Times=AP_Param(:,3);
                AP_Peak_Vm=AP_Param(:,4);
                AP_Duration=AP_Param(:,6);
                
                [Vm_Sub]=Function_CutAPs(MembranePotential, SR_Vm, AP_Peak_Times, AP_Thrs_Times, AP_Thrs_Vm, AP_Peak_Vm);
                
            else
                Vm_Sub=MembranePotential;
            end
            
            % update 1 Cell variables from the current sweep
            
            Tot_AP_Thrs=vertcat(Tot_AP_Thrs, AP_Thrs_Vm);
            Tot_AP_Dur=vertcat(Tot_AP_Dur, AP_Duration);
            
            Mean_Vm=[];
            SD_Vm=[];
            [Mean_Vm, SD_Vm]=Function_SubThrsVm(Vm_Sub, SR_Vm, TimeWindow); % Compute mean Vm and Vn SD
            
            Tot_Mean_Vm=vertcat(Tot_Mean_Vm, Mean_Vm);
            Tot_SD_Vm=vertcat(Tot_SD_Vm, SD_Vm);
            
            FFT_Mtrx=[];
            [FFT_Mtrx]=Function_Compute_FFTs(Vm_Sub, SR_Vm, TimeWindow); % Compute Vm FFT
            Tot_FFT=horzcat(Tot_FFT, FFT_Mtrx);
            
            
        end
        
        % Fill result structure with results
        
        result.(cell2mat(Cell_Types(tp))).Firing_Rate(c,1)=Tot_Numb_APs/Tot_Rec_Dur;
        result.(cell2mat(Cell_Types(tp))).AP_Threshold(c,1)=mean(Tot_AP_Thrs, 'omitnan');
        result.(cell2mat(Cell_Types(tp))).AP_Duration(c,1)=mean(Tot_AP_Dur, 'omitnan');
        result.(cell2mat(Cell_Types(tp))).mean_Vm(c,1)=mean(Tot_Mean_Vm);
        result.(cell2mat(Cell_Types(tp))).SD_Vm(c,1)=mean(Tot_SD_Vm);
        
        mean_FFT=mean(Tot_FFT,2);
        result.(cell2mat(Cell_Types(tp))).FFT{c,:}=mean_FFT;
             
        Step=TimeWindow*SR_Vm;
        nfft = 2^nextpow2(Step);
        FreqBandppt=round((FreqBandLim*nfft)/SR_Vm)+1;
        
        result.(cell2mat(Cell_Types(tp))).FFT_LF(c,1)=mean(mean_FFT(FreqBandppt(1):FreqBandppt(2)));
        result.(cell2mat(Cell_Types(tp))).FFT_HF(c,1)=mean(mean_FFT(FreqBandppt(3):FreqBandppt(4)));
                
    end
    
    
end

%% SAVE THE RESULT STRUCTURE

disp('SAVING RESULTS')

StructureName='Result_1';
save([PathSaveResults filesep StructureName], 'result','-v7.3');

disp('RESULT SAVED')