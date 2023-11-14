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

Pre_Window=0.5; % time before whisking onset (s)
Post_Window=0.5; % time after whisking onset (s)
Min_Event_Dur=0.2; % minimal duration of whisking episode to be considered
Min_ITI=0.5; % minimal interval between 2 consecutive whisking episode (s)
Min_Numb_Trial=3; % minimal number of trial
Vm_Levels_Times=[-0.5 -0.3; 0 0.2]; % Times to compute changes in Vm
AP_FR_Times=[-0.5 -0.3; 0 0.2]; % Times to compute changes in Firing Rate


%% Select for one Sweep type
dataSelect=[];
myFields = fieldnames(Data);
for thisField = 1:length(myFields)
    dataSelect.(cell2mat(myFields(thisField))) = Data.(cell2mat(myFields(thisField)))(strcmp(Data.Sweep_Type, Sweep_Type),:);
end

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
        
        result.(cell2mat(Cell_Types(tp))).Cell_Name{c,1}=data1Cell.Cell_ID{1,1};
        result.(cell2mat(Cell_Types(tp))).Cell_Type{c,1}=data1Cell.Cell_Type{1,1};
        result.(cell2mat(Cell_Types(tp))).Cell_Depth(c,1)=data1Cell.Cell_Depth(1,1);
        result.(cell2mat(Cell_Types(tp))).Cell_Layer{c,1}=data1Cell.Cell_Layer{1,1};
        result.(cell2mat(Cell_Types(tp))).Cell_Anatomy{c,1}=data1Cell.Cell_Anatomy{1,1};
        result.(cell2mat(Cell_Types(tp))).Cell_Anatomy{c,2}=data1Cell.Cell_Anatomy{1,2};
        result.(cell2mat(Cell_Types(tp))).SR_Vm(c,1)=data1Cell.Sweep_MembranePotential_SamplingRate(1,1);
        result.(cell2mat(Cell_Types(tp))).SR_WP(c,1)=data1Cell.Sweep_WhiskerAngle_SamplingRate(1,1);
        
        Tot_Vm_Mtrx=[];
        Tot_APs_Mtrx=[];
        Tot_WP_Mtrx=[];
        Event_Times=[];
        
        for sweep=1:size(data1Cell.Sweep_Counter,1) % loop through the sweeps from 1 cell
            
            MembranePotential=[];
            SR_Vm=[];
            Rec_Dur=[];
            
            Event_Times=data1Cell.Sweep_WhiskingTimes{sweep,1};
  
            MembranePotential=data1Cell.Sweep_MembranePotential{sweep,1};
            SR_Vm=data1Cell.Sweep_MembranePotential_SamplingRate(sweep,1);
            
            WhiskerAngle=data1Cell.Sweep_WhiskerAngle{sweep,1};
            SR_WP=data1Cell.Sweep_WhiskerAngle_SamplingRate(sweep,1);
            
            AP_Vm_Deriv_Thrs=data1Cell.Cell_APThreshold_Slope(1,1);
            [AP_Param]=Function_Detect_APs(MembranePotential, SR_Vm, AP_Vm_Deriv_Thrs);
            
            AP_Peak_Times=[];
            AP_Peak_Vm=[];
            AP_Thrs_Times=[];
            AP_Thrs_Vm=[];
            APs_Mtrx=[];
             
            Vm_Sub=[];
            AP_vect=[];
            
            if ~isempty(AP_Param)
                
                AP_Peak_Times=AP_Param(:,3);
                AP_Thrs_Times=AP_Param(:,1);
                AP_Thrs_Vm=AP_Param(:,2);
                AP_Peak_Vm=AP_Param(:,4);
               
                
                [Vm_Sub]=Function_CutAPs(MembranePotential, SR_Vm, AP_Peak_Times, AP_Thrs_Times, AP_Thrs_Vm, AP_Peak_Vm);
                
                Length_Vect=length(Vm_Sub);
                [AP_vect]=Function_Times2Vect(AP_Peak_Times, SR_Vm, Length_Vect);
            else
                Vm_Sub=MembranePotential;
                AP_vect(1:size(Vm_Sub,1),1)=0;
            end
     
            Vm_Mtrx=[];
            [Vm_Mtrx]=Function_Event_Triggered_Signal(Vm_Sub, SR_Vm, Event_Times, Pre_Window, Post_Window, Min_Event_Dur, Min_ITI);
            Tot_Vm_Mtrx=horzcat(Tot_Vm_Mtrx, Vm_Mtrx);
            
            WP_Mtrx=[];
            [WP_Mtrx]=Function_Event_Triggered_Signal(WhiskerAngle, SR_WP, Event_Times, Pre_Window, Post_Window, Min_Event_Dur, Min_ITI);
            Tot_WP_Mtrx=horzcat(Tot_WP_Mtrx, WP_Mtrx);
            
            APs_Mtrx=[];
            [APs_Mtrx]=Function_Event_Triggered_Signal(AP_vect, SR_Vm, Event_Times, Pre_Window, Post_Window, Min_Event_Dur, Min_ITI);
            Tot_APs_Mtrx=horzcat(Tot_APs_Mtrx, APs_Mtrx);
           
          
        end
 %%       
        if size(Tot_Vm_Mtrx,2)>=Min_Numb_Trial;
            Vm_avg=[];
            WP_avg=[];
            AP_avg=[];
            Vm_avg=mean(Tot_Vm_Mtrx,2)';
            WP_avg=mean(Tot_WP_Mtrx,2)';
            AP_avg=mean(Tot_APs_Mtrx,2)';
            
            result.(cell2mat(Cell_Types(tp))).Vm_avg{c,1}=Vm_avg;
            result.(cell2mat(Cell_Types(tp))).WP_avg{c,1}=WP_avg;
            result.(cell2mat(Cell_Types(tp))).AP_avg{c,1}=AP_avg;
            result.(cell2mat(Cell_Types(tp))).Numb_Onset(c,1)=size(Tot_Vm_Mtrx,2);
            
            for t=1:size(Vm_Levels_Times,1)
                Vm_lev1=round((Vm_Levels_Times(t,1)+Pre_Window)*SR_Vm+1);
                Vm_lev2=round((Vm_Levels_Times(t,2)+Pre_Window)*SR_Vm);
                
                result.(cell2mat(Cell_Types(tp))).Vm_Amplitude(c,t)=mean(Vm_avg(1,Vm_lev1:Vm_lev2));
                
            end
            
            for t=1:size(AP_FR_Times,1)
                AP_lev1=round((AP_FR_Times(t,1)+Pre_Window)*SR_Vm+1);
                AP_lev2=round((AP_FR_Times(t,2)+Pre_Window)*SR_Vm);
                
                result.(cell2mat(Cell_Types(tp))).AP_FiringRate(c,t)=sum(AP_avg(1,AP_lev1:AP_lev2))/(AP_FR_Times(t,2)-AP_FR_Times(t,1));
                
            end
            
            
        else
            result.(cell2mat(Cell_Types(tp))).Vm_avg{c,1}=NaN;
            result.(cell2mat(Cell_Types(tp))).WP_avg{c,1}=NaN;
            result.(cell2mat(Cell_Types(tp))).AP_avg{c,1}=NaN;
            result.(cell2mat(Cell_Types(tp))).Numb_Onset(c,1)=size(Tot_Vm_Mtrx,2);
            result.(cell2mat(Cell_Types(tp))).Vm_Amplitude(c,1:size(Vm_Levels_Times,1))=NaN;
            result.(cell2mat(Cell_Types(tp))).AP_FiringRate(c,1:size(AP_FR_Times,1))=NaN;
            
        end
        
        
    end
    
    
end

%% SAVE THE RESULT STRUCTURE

disp('SAVING RESULTS')

StructureName='Result_2';
save([PathSaveResults filesep StructureName], 'result','-v7.3');

disp('RESULT SAVED')