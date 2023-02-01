function [Vm_Sub]=Function_CutAPs(Vm, SR_Vm, AP_Peak_Times, AP_Thrs_Times, AP_Thrs_Vm, AP_Peak_Vm)

% This function computes the subthreshold Vm ('Vm_Sub') by cutting the APs from the raw Vm signal

%% INPUTS:
% Vm = membrane potential vector (V)
% SR_Vm = sampling rate of the Vm vector (sample/s)
% AP_Peak_Times = vector containing the time (s) of the peak of each AP
% AP_Thrs_Times = vector containing the time (s) of the threshold of each AP
% AP_Thrs_Vm = vector containing the Vm (V) of the peak of each AP
% AP_Peak_Vm = vector containing the Vm (V) of the threshold of each AP

%% OUPUT:
% Vm_Sub= subthreshold Vm vector (V)


%% Parameters

AP_End_Wind=0.015; % time window to look for AP end after AP peak (s)

%%

Vm_Sub=[];
Vm_Sub=Vm;

AP_Thrs_Index=round(AP_Thrs_Times.*SR_Vm);
AP_Peak_Index=round(AP_Peak_Times.*SR_Vm);

for Ind=1:size(AP_Peak_Index,1)
    
    pt1=AP_Thrs_Index(Ind,1);
    pt2=AP_Peak_Index(Ind,1);
    pt4=min(length(Vm),pt2+round(AP_End_Wind*SR_Vm));
    
    Vm_Thrs=Vm(pt1,1);
    
    AP_seg=smooth(Vm(pt1:pt4,1),SR_Vm/2000);
    
    cond=0;
    n2=pt2-pt1+20;
    
    while cond==0
        
        n2=n2+1;
        
        if n2>length(AP_seg)
            pt_end=pt1+n2;
            cond=1;
        else
            D_Vm1=AP_seg(n2-1,1)-Vm_Thrs;
            D_Vm2=AP_seg(n2,1)-Vm_Thrs;
            
            if D_Vm2<0
                pt_end=pt1+n2;
                cond=1;
            end
            
            if D_Vm2>D_Vm1
                pt_end=pt1+n2;
                cond=1;
            end
        end
        
    end
    
    pt3=min(length(Vm), pt_end);
    
    % make a segment 'in' between Vm(pt1) and Vm(pt2)
    Delta_Vm=Vm(pt3)-Vm(pt1);
    in=0:1:pt3-pt1; % create a small vector starting at 0, with increment of 1, of (pt2-pt1) points
    in=(in./(pt3-pt1))*Delta_Vm; % create a segment of (pt2-pt1) points from 0 to Delta_Vm
    in=in+Vm(pt1); % add the Vm at pt1 to the segment in
    
    Vm_Sub(pt1:pt3,1)=in;  
    
end


end