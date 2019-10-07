function implay_RR( DATA, colmap, limits)
%Coloured and scaled implay: limits = [min max], colmap = {jet, gray..etc}

handle = implay(DATA);

if nargin == 1
    colmap = 'gray';
    limits(1) = min(min(min(DATA)));
    limits(2) = 0.75* max(max(max(DATA)));
elseif nargin == 2
    limits(1) = min(min(min(DATA)));
    limits(2) = 0.75* max(max(max(DATA)));
end

handle.Visual.ColorMap.UserRangeMin = limits(1);
handle.Visual.ColorMap.UserRangeMax = limits(2);
handle.Visual.ColorMap.UserRange = 1;

c = [colmap '(256)'];
handle.Visual.ColorMap.MapExpression = c;

set(findall(0,'tag','spcui_scope_framework'),'position',[1080/2 -150 900 900]);

end

