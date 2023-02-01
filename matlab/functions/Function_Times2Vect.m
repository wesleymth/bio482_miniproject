function [AP_vect]=Function_Times2Vect(AP_times, SR_Vm, Length_Vect);

% This function generate a binary vector of of the APs with same resolution
% and length as the original Membrane potential signal

%% INPUTS:
% AP_times = vector containing the times of each AP (s)
% SR_Vm = sampling rate of the original Vm trace (sample / s)
% Length_Vect = number of point of the original Vm trace

%% OUTPUT:
% AP_vect = binary vector (0 = non AP, 1= AP) with resolution 'SR_Vm' and
% length 'Length_Vect'

%%
AP_vect=[];

AP_vect(1:Length_Vect,1)=0;

if ~isnan(AP_times)
    
    for t=1:size(AP_times,1)
        
        AP_pt=round(AP_times(t,1)*SR_Vm);
        AP_vect(AP_pt,1)=1;
        
    end
end


end