function plotVertLine(values,modifiers)
% Plots a vertical X=* line at each of the entered values

V = axis;

currHold = ishold;
hold on

for i = 1:length(values)
    if ~exist('modifiers','var')
        plot([values(i) values(i)],[V(3) V(4)])
    else
        plot([values(i) values(i)],[V(3) V(4)],modifiers{:})
    end
    
end

axis(V)

if ~currHold
    hold off
end