function [FFT_Mtrx]=Function_Compute_FFTs(Vm_Sub, SR_Vm, TimeWindow)

% this function generate a matrix containing the FFTs computed for segments
% of duration 'TimeWindow' cut from the subthreshold Vm (Vm_Sub)

%% INPUTS:
% Vm_Sub = subthreshold membrane potential vector (V)
% SR_Vm = sampling rate of the Vm vector (sample/s)
% TimeWindow = duration of the time window used to cut Vm (s)

%% OUPUT:
% FFT_Mtrx= Matrix containing the FFT computed for each time window

%%

Numb_Wind=floor((length(Vm_Sub)/SR_Vm)/TimeWindow);

step=TimeWindow*SR_Vm;
nfft = 2^nextpow2(step); % Calculate the next power of 2 for padding
window = hanning(step); % Hanning Window

for s=1:Numb_Wind
    pt1=step*(s-1)+1;
    pt2=step*s;
    Seg=Vm_Sub(pt1:pt2);
    Seg=Seg-mean(Seg); % remove the DC component at 0 Hz
    Seg_fft=fft(Seg.*window,nfft); % compute the FFT for the segment Seg
    
    P2=abs(Seg_fft/step); % compute the spectrogram
    P1=P2(1:nfft/2+1);
    P1(2:end-1)=2*P1(2:end-1);
    FFT_Mtrx(:,s)=P1; % make a matrix containing all the FFTs
    
end

end




        
 