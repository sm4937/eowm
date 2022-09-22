function [cw_biases,meridian_biases] = get_angular_biases(angles_A,angles_B)
%get_angular_distance This functiion pulls how biased people are
%   towards counter- or clockwise directions, or towards the horizontal
%   or vertical meridians.
% Angles_A should be the subject's estimate or response
% Angles_B should be the true answer
% Either Angles_A or Angles_B can be a vector or a single number. 

% transformed_angs = 1:360; transformed_angs = transformed_angs - [zeros(1,180) 360*ones(1,180)];
% angles_A = transformed_angs(angles_A+1); % index in to this here using original angles_A vector
% angles_B = transformed_angs(angles_B+1);

% angles come in from 0 to 359, transform to 1 to 360, then 
% 1 to 180, -179 to 0

% if row vectors, make column vectors
if size(angles_A,2)>size(angles_A,1)
    angles_A = angles_A';
end
if size(angles_B,2)>size(angles_B,1)
    angles_B = angles_B';
end

raw_distance = angles_A - angles_B;
[angular_distance,which] = min([abs(raw_distance) abs(360-raw_distance)],[],2);

% clockwise is 2
% counter-clockwise is 1
% just like the button-presses

cw_biases = 2.*ones(length(angles_A),1);
cw_biases(which==2,:) = 1;
cw_biases(raw_distance<0,:) = 1;

meridian_biases = NaN(length(angles_A),4);
% column 1: degrees of bias towards 0 (360) degrees
% column 2: degrees of bias towards 90 degrees
% column 3: degrees of bias towards 180 degrees
% column 4: degrees of bias towards 270 degrees

% first: identifty quadrant of true stimulus
quadrants = ceil((angles_B+1)./90);
%raw_distance will be negative if the true answer is BIGGER
%than the decoded answer (closer to cw meridian)
%it will be positive if true answer is smaller (closer to ccw meridian)
for err = 1:length(angles_B)
    
    
    real_angle = angles_B(err);
    decoded_angle = angles_A(err);
    
    real_distances = abs(real_angle - [360 90 180 270]);
    real_distances(1) = min(abs([real_angle - 0, real_angle-360]));
    decoded_distances = abs(decoded_angle - [360 90 180 270]);
    decoded_distances(1) = min(abs([decoded_angle - 0, decoded_angle-360]));

    % if the decoded distance is SMALLER than the real distance
    % that's a positive bias TOWARDS that meridian
    
    meridian_biases(err,:) = real_distances - decoded_distances;
    %so positive numbers here are drift TOWARDS
    
    
end



end

