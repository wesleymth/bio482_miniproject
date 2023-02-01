function Function_Plot_Mean_FFT(FFT_Matrix, SR, color)

% This function compute and plot the mean and sem of the FFT matrix on a
% semi-log scale

%% INPUTS:
% FFT_Matrix = matrix containing the individual FFTs
% SR = sampling rate of the original signal (sample / s)
% color = RGB color vector [R G B]

%%


FFT_Mean=mean(FFT_Matrix,1);
FFT_sem=std(FFT_Matrix,1)/sqrt(size(FFT_Matrix,1));

FFT_Mean_sup=FFT_Mean+FFT_sem;
FFT_Mean_inf=FFT_Mean-FFT_sem;

f=[];
nfft = (size(FFT_Matrix,2)-1)*2;
f = SR*(0:(nfft/2))/nfft;

semilogx(f,FFT_Mean, 'Color', color, 'Linewidth', 2)
hold on
semilogx(f,FFT_Mean_sup, 'Color', color, 'Linewidth', 0.5)
hold on
semilogx(f,FFT_Mean_inf, 'Color', color, 'Linewidth', 0.5)
hold on

end