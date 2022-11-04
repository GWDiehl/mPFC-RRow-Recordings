% Plots lines @ X=0 and Y=0 on the current plot to provide quadrant
% designations
% Input finish identifies the end point of each set of lines


function plotQuadrantBounds(modifiers)

V = axis;
currHold = ishold;

dimentions = length(V)/2;

hold on
if ~exist('modifiers','var')
    if dimentions == 2
        plot([V(1) V(2)],[0 0])
        plot([0 0],[V(3) V(4)])
    elseif dimentions == 3
        plot3([V(1) V(2)],[0 0],[0 0])
        plot3([0 0],[V(3) V(4)],[0 0])
        plot3([0 0],[0 0],[V(5) V(6)])
    end
else
    if dimentions == 2
        plot([V(1) V(2)],[0 0],modifiers{:})
        plot([0 0],[V(3) V(4)],modifiers{:})
    elseif dimentions == 3
        plot3([V(1) V(2)],[0 0],[0 0],modifiers{:})
        plot3([0 0],[V(3) V(4)],[0 0],modifiers{:})
        plot3([0 0],[0 0],[V(5) V(6)],modifiers{:})
    end
end

axis(V)

if ~currHold
    hold off
end
