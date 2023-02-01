function [AP_PSTH]=function_PSTH(Onset_AVG_AP,SR_Vm, Pre_Window, Post_Window, bin_size)

% This function generates a peri-stimulus time histogram (PSTH) based on the averaged AP signal around event time 

%% INPUTS:
% Onset_AVG_AP = vector of the averaged AP signal around event onset time
% SR_Vm = Sampling rate of the membrane potential (sample / s)
% Pre_Window = time before event onset time (s)
% Post_Window = time after event onset time (s)
% bin_size = size of the bins to compute the PSTH (s)

%% OUPUT:
% AP_PSTH = matrix containing the time vector (column 1) (s) and the firing
% rate (column 2) (Hz) for each time bin.

%% 
AP_PSTH=[];
AP_PSTH_pre=[];
AP_PSTH_post=[];

bin_pt=round(bin_size*SR_Vm);

pt0=floor(-1*Pre_Window*SR_Vm);
pt2=pt0;
pt_min=pt0-floor(pt0/bin_pt)*bin_pt+1;

cnt=0;

while(pt2>pt_min)
    
    pt1=pt0-(bin_pt*(cnt));
    pt2=pt1-bin_pt+1;
    
    AP_PSTH_pre(cnt+1,1)=-1*bin_size*(cnt+1);
    AP_PSTH_pre(cnt+1,2)=sum(Onset_AVG_AP(pt2:pt1,1))/bin_size;
    cnt=cnt+1;
    
end

pt0=floor(-1*Pre_Window*SR_Vm);
pt2=pt0;
pt_max=pt0+floor((Post_Window*SR_Vm)/bin_pt)*bin_pt;

cnt=0;

while(pt2<pt_max)
    
    pt1=pt0+(bin_pt*cnt);
    pt2=pt1+bin_pt;
    
    AP_PSTH_post(cnt+1,1)=bin_size*(cnt);
    AP_PSTH_post(cnt+1,2)=sum(Onset_AVG_AP(pt1:pt2,1))/bin_size;
    
    cnt=cnt+1;
    
end

AP_PSTH=vertcat(AP_PSTH_pre, AP_PSTH_post);

AP_PSTH=sortrows(AP_PSTH);

end