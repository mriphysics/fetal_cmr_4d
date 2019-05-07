function [c,ceq] = nonlconfn(x)
%NONLCONFN  Nonlinear constraints function.
%
%   [c,ceq] = NONLCONFN(x) returns nonlinear inequality constraints c and 
%   nonlinear equality constraints ceq of the form c(x) <= 0 and ceq(x) = 0.
%
%   See also CARDSYNC_INTERSLICE.

%   jfpva (joshua.vanamerom@kcl.ac.uk)

c   = zeros(size(x));
ceq = zeros(size(x));

end   % nonlconfn(...)