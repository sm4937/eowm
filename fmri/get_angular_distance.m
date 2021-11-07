function [angular_distance] = get_angular_distance(angles_A,angles_B)
%get_angular_distance This functiion pulls the true circular distance of
%pairs of angles
%   angles_A and angles_B can be lists of angles, as long as they're the
%   same size, or just one number each 

raw_distance = abs(angles_A - angles_B);
angular_distance = min([raw_distance abs(360-raw_distance)],[],2);


end

