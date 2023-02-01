function [AP_Param]=Function_Detect_APs(MembranePotential, SR_Vm, Vm_Deriv_Thrs)

% This function detect the APs from the membrane potential signal and
% compute different parameters for each AP. The results are saves in the
% matrix 'AP_Param'.

%% INPUTS:
% MembranePotential = vector containing Vm data (V)
% SR_Vm = sampling rate of the MembranePotential (sample / s)
% Vm_Deriv_Thrs = threshold to detect AP initiation (V/s)

%% OUPUT:
% AP_Param = a matrix that contains the measured parameters for each detected AP
% column 1 = AP Threshold time (s)
% column 2 = AP Threshold Vm (V)
% column 3 = AP peak time (s)
% column 4 = AP peak Vm (V)
% column 5 = AP amplitude (V)
% column 6 = AP duration at 1/2 amplitude (ms)


%% Parameters

AP_Win=0.0015; % time (s) to search for AP's peak
Min_AP_Amp= 0.005; % Minimum AP amplitude (V)
AP_length=round(AP_Win*SR_Vm);

%%

AP_Thrs_Index=[];
AP_Param=[]; % 1col= Thrs_Times, 2col=Thrs_Vm, 3col=Peak_Times, 4col=Peak_Vm, 5col= Peak_Amp, 6col=AP_Dur

Vm_Deriv=diff(MembranePotential)*SR_Vm;
AP_Thrs_Onset= diff((Vm_Deriv-Vm_Deriv_Thrs)./abs(Vm_Deriv-Vm_Deriv_Thrs));

Vm_Med=median(MembranePotential);

[Peaks, AP_Thrs_Index] = findpeaks(AP_Thrs_Onset, 'MinPeakHeight', 0.1, 'MinPeakProminence',0.5, 'MinPeakDistance',SR_Vm*0.001);

if isempty(AP_Thrs_Index)==0
    
    AP_cnt=1;
    
    for i=1:size(AP_Thrs_Index,1)
        
        pt1=AP_Thrs_Index(i,1);
        pt2=AP_Thrs_Index(i,1)+AP_length;
        
        if pt2<length(MembranePotential)
            
            AP_Seg=MembranePotential(pt1:pt2,1);
            [Max, Ind]=max(AP_Seg);
            
            if ~isempty(Ind)
                
                AP_Index=pt1+Ind(1,1)-1;
                
                AP_Amp=MembranePotential(AP_Index,1)-MembranePotential(pt1,1);
                
                if AP_Amp>Min_AP_Amp && MembranePotential(AP_Index,1)>Vm_Med
                    
                    AP_Param(AP_cnt,1)=AP_Thrs_Index(i,1)/SR_Vm; % Thrs Time
                    AP_Param(AP_cnt,2)=MembranePotential(AP_Thrs_Index(i,1),1); % Thrs Vm
                    AP_Param(AP_cnt,3)=AP_Index/SR_Vm; % Peak Time
                    AP_Param(AP_cnt,4)=MembranePotential(AP_Index,1); % Peak Vm
                    AP_Param(AP_cnt,5)=AP_Amp; % Peak Amp

                    Ind=AP_Thrs_Index(i,1);
                    pt2=AP_Index;
                    % we define a time window for the AP ...
                    pt1=pt2-0.002*SR_Vm; % ... 2 ms before the peak ...
                    pt3=pt2+0.003*SR_Vm; % ... and 3 ms after the peak.
                    
                    if pt1>0 && pt3<length(MembranePotential)
                                              
                        Vm_HalfAmp=MembranePotential(Ind,1)+AP_Amp/2; % Vm at half amplitude
                        
                        sAP_Seg=MembranePotential(pt1:pt3,1); % cut a segment of the VM that contains the AP
                        sAP_Seg=sAP_Seg-Vm_HalfAmp; % substract the Vm at half-amplitude
                        sAP_OnOff=diff(sAP_Seg./abs(sAP_Seg)); % compute the binary signal
                        
                        [sAP_Max, sAP_Indmax]=max(sAP_OnOff); % identify index begening AP at half amplitude
                        [sAP_Min, sAP_Indmin]=min(sAP_OnOff); % identify index end AP at half amplitude
                        
                        AP_Param(AP_cnt,6)=((sAP_Indmin-sAP_Indmax)/SR_Vm)*1000; % compute duration at half-amplitude
                        
                        
                    else
                        AP_Param(AP_cnt,6)=NaN;
                    end
                    
                    AP_cnt=AP_cnt+1;
                    
                end
            end
        end
    end
    
end

% Remove outliers = APs that have peak height and Amplitude < 5x std

if ~isempty(AP_Param)
    
    Amp_Lim_Inf=min(median(AP_Param(:,5))-5*std(AP_Param(:,5)), 0.03);
    Peak_Lim_Inf= min(median(AP_Param(:,4))-5*std(AP_Param(:,4)), -0.02);
    
    
    cnt_Max=size(AP_Param,1);
    
    for i=1:cnt_Max
        
        cnt=cnt_Max-i+1;
        if AP_Param(cnt,5)< Amp_Lim_Inf && AP_Param(cnt,4)< Peak_Lim_Inf
            AP_Param(cnt,:)=[];
        end
    end
    
end

end