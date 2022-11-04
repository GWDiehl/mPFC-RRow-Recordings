

function plotIDLine(modifiers)
% Plots an Identity x=y plane on the current plot

% GWDiehl

V = axis;
currHold = ishold;

pt1 = min(V(1:2:end));
pt2 = max(V(2:2:end));

dimentions = length(V)/2;

hold on
if ~exist('modifiers','var')    
    if dimentions == 2
        plot([pt1 pt2], [pt1 pt2])
    elseif dimentions == 3
        plot3([pt1 pt2],[pt1 pt2],[pt1 pt2])
    end
else
     if dimentions == 2
        plot([pt1 pt2], [pt1 pt2],modifiers{:})
    elseif dimentions == 3
        plot3([pt1 pt2],[pt1 pt2],[pt1 pt2],modifiers{:})
     end
end

axis(V)

if ~currHold
    hold off
end