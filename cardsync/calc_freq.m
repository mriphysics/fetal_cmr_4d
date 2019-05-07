function f = calc_freq( nFrame, dt )
%CALC_FREQ  Calculate x-f frequency.
% 
%   f = CALC_FREQ( nFrame, dt ) returns a vector of frequencies for x-f
%   space given number of frames, nFrame, and frame duration, dt, in
%   seconds.

%   jfpva (joshua.vanamerom@kcl.ac.uk) 


iFdc = floor(nFrame/2)+1;

fMax = 1/dt;

f = fMax * ( (1:nFrame) - iFdc ) / nFrame;


end  % calc_freq(...)

