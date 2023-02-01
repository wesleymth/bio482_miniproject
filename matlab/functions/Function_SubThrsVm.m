
function [Mean_Vm, SD_Vm]=Function_SubThrsVm(MembranePotential, SR_Vm, TimeWindow)

% This function computes the mean and SD of the Vm for consecutive time
% windows

%% INPUTS = 
% MembranePotential = vector containing Vm data (V)
% SR_Vm = sampling rate of the MembranePotential (sample / s)
% TimeWindow = duration of the time window used to cut Vm (s)

%% OUPUTS =
% Mean_Vm = vector containing the mean Vm computed for each time window (V)
% SD_Vm = vector containing the SD of the Vm computed for each time window (V)

%%

Numb_Wind=floor((length(MembranePotential)/SR_Vm)/TimeWindow);

for window=1:Numb_Wind
    
    pt1=1+TimeWindow*SR_Vm*(window-1);
    pt2=pt1+TimeWindow*SR_Vm-1;
    
    Mean_Vm(window,1)=mean(MembranePotential(pt1:pt2,1));
    SD_Vm(window,1)=std(MembranePotential(pt1:pt2,1));
        
end


end
