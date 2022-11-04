function DesiredLocation = rescaleIn3D(DesiredDist,TotalDist,Start,End)

% Takes a Start and End location in arbitrary 3D space (or ND really),
% scales such that the the distance between the points is TotalDistance,
% and localizes a point along the vector that is DesiredDistance.

% Note that this takes the original start location and re-computes the desired
% END location along the original vector


% Originally written for localizing unit data in standardized 3D space but
% could be used for anything.

% GWD March 2020

origOffset = Start;
ShiftedEnd = End - origOffset; % Move so it starts at the origin

% What proportion of the total distance are we going
scaleFactor = DesiredDist/TotalDist;


ScaledEnd = ShiftedEnd*scaleFactor; % Make it the right length

DesiredLocation = ScaledEnd + origOffset; % Put it back to the right start location




