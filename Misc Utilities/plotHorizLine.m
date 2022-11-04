function plotHorizLine(values,modifiers)
% Plots a horizontal Y=* line at each of the entered values

V = axis;

currHold = ishold;    
    hold on
for i = 1:length(values)
    
    if ~exist('modifiers','var')
        plot([V(1) V(2)],[values(i) values(i)])
    else
        plot([V(1) V(2)],[values(i) values(i)],modifiers{:})
    end
    
end

axis(V)
if ~currHold
    hold off
end