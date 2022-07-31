%Size of training band in range dimension = 8
%Size of guard band in range dimension = 4
%Size of training band in doppler dimension = 4
%Size of guard band in the doppler dimension = 2

result = (12*2+1) * (6*2+1) - (4*2+1)*(2*2+1);

%Steering Angle
%Determine the steering angle of an antenna beam in degrees for the given phase increment:
% Frequency of operation = 77 GHz 
% Speed of light = 3*10^8 m/s 
% Phase increment = 45 degrees 
% Antenna element spacing = wavelength/2

f = 77e9;
c = 3e8;
lambda = c/f;
d = lambda / 2;
% pi / 4 = 360 * d * sin (theta) / lambda;
theta = asin (lambda * (45) / (360 * d));
theta = theta * 180 / pi