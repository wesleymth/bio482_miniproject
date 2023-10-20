function [Signal_Mtrx]=Function_Event_Triggered_Signal(Signal, SR, Event_Times, Pre_Window, Post_Window, Min_Event_Dur, Min_ITI)

% This function generate a matrix of segments cut from the 'Signal' around ('Pre_Window' before and 'Post_Window' after)
% each time event 

%% INPUTS:
% Signal = signal to be cut into segments around event times
% SR = sampling rate (sample /s)
% Event_Times = 2xN Matrix containing the onset and offset times of each event (s)
% Pre_Window = time before event time to cut signal (s)
% Post_Window = time after event time to cut signal (s)
% Min_Event_Dur = minimal duration of the event to consider (s)
% Min_ITI = minimal interval between the onset of the event and the offset of the preceding event (s)

%% OUPUT:
% Signal_Mtrx = Matrix containing the segments cut around event times

%%

Signal_Mtrx=[];
cnt=1;

for i=1:size(Event_Times, 1)
    
    Event_Dur=[];
    ITI=[];
    % compute the event duration
    Event_Dur=Event_Times(i,2)-Event_Times(i,1);
    % compute the ITI (ITI = event time for the 1st event)
    if i==1
        ITI=Event_Times(i,1);
    else
        ITI=Event_Times(i,1)- Event_Times(i-1,2);
    end
    
    if Event_Dur> Min_Event_Dur && ITI>Min_ITI
       
        pt1=ceil((Event_Times(i,1)-Pre_Window)*SR); 
        pt2=pt1+ceil((Pre_Window+Post_Window)*SR)-1;
        
        if pt1>0 && pt2<length(Signal)
            
            Signal_Mtrx(:,cnt)=Signal(pt1:pt2,1); % cut the Vm around the event time
            
            cnt=cnt+1;
        end
    end
    
end



end