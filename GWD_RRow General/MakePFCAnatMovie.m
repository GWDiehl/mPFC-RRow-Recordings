
function h = MakePFCAnatMovie(outputDir,fileName,varargin)

anatGroup = 'No_VO';

process_varargin(varargin);

fullFN = [outputDir,fileName];

[~, ~, ~, ~, ~,h] = LocalizeUnits_mPFC('anatGroup',anatGroup,'anatplotType','Both');

axis off
view(0,0)
axis tight manual % this ensures that getframe() returns a consistent size


stepSize = 1;
pauseTime = .025;
azSteps = stepSize:stepSize:360;
nSteps = length(azSteps);
elSteps = [linspace(0,35,nSteps/4) linspace(35,0,nSteps/4)];
elSteps = repmat(elSteps,1,2);


pauseDur = 10;

% Capture the plot as an image
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,fullFN,'gif', 'Loopcount',inf);


% Capture Paused Frame
for n = 1:pauseDur
%     pause(pauseTime)
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    imwrite(imind,cm,fullFN,'gif','WriteMode','append');
end

% XAxis Rotation
for n = 1:nSteps
view(azSteps(n),0)
% pause(pauseTime)

    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    imwrite(imind,cm,fullFN,'gif','WriteMode','append');
end

% Capture Paused Frame
for n = 1:pauseDur
%     pause(pauseTime)
    
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    imwrite(imind,cm,fullFN,'gif','WriteMode','append');
end

% Oblique Rotation
for n = 1:nSteps
view(azSteps(n),elSteps(n))
% pause(pauseTime)

    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    imwrite(imind,cm,fullFN,'gif','WriteMode','append');
end


