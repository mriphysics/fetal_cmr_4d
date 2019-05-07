function [ timeSinceTrigger, cardiacPhase, rrInterval ] = calc_cardiac_timing( tRealtime, tTrigger )
%CALC_CARDIAC_TIMING  Calculate cardiac cine timing.
% 
%   CALC_CARDIAC_TIMING( tRealtime, tTrigger ) calculates cardiac cine 
%   timing given realtime image sequence timing, tRealtime, and cardiac 
%   trigger times, tTrigger.
% 
%   [timeSinceTrigger,cardiacPhase,rrInterval] = CALC_CARDIAC_TIMING(...)
%   returns time since last trigger, cardiac phase, and R-R interval for
%   each frame in the realtime image sequence.

%   jfpva (joshua.vanamerom@kcl.ac.uk)


%% Check Input

assert( min(tRealtime)>=min(tTrigger),...
    'need trigger time before first realtime image frame'); 


assert( max(tRealtime)<=max(tTrigger),...
    'need trigger time after last realtime image frame'); 


%% Setup

nFrame           = length( tRealtime );
cardiacPhase     = nan( size( tRealtime ) );
timeSinceTrigger = nan( size( tRealtime ) );
rrInterval       = nan( size( tRealtime ) );


%% Calculate Cardiac Timing

for iFrame = 1:nFrame
    
    tNextTrigger = min( tTrigger( ( tTrigger >  tRealtime(iFrame))));
    tPrevTrigger = max( tTrigger( ( tTrigger <= tRealtime(iFrame))));
    
    timeSinceTrigger(iFrame) = tRealtime(iFrame) - tPrevTrigger;
    
    rrInterval(iFrame)       = tNextTrigger - tPrevTrigger;
    
    cardiacPhase(iFrame)     = timeSinceTrigger(iFrame) / rrInterval(iFrame);

end


end  % calc_cardiac_timing(...) 

