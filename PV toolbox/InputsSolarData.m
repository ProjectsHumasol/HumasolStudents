% This file store all the required user defined parameters to run the
% "DetermineSolarData" script

% Tilt angle of solar panels (generally constant for a given country, in
% degrees)
panelTiltAngle = 6;

% Solcast Data
file = "Data_Solcast_PT15M15-16.xlsx";
temperatureRange = "D97:D70272";
azimuthRange = "E97:E70272";
DHIRange = "H97:H70272";
DNIRange = "I97:I70272";
GHIRange = "K97:K70272";
zenithRange = "R97:R70272";
[DHI, DNI, GHI, zenith, azimuth,temperature] = readSolcastData(file, DHIRange, DNIRange, GHIRange, zenithRange, azimuthRange, temperatureRange);
PTM = 4; % Number of measurements per hour, determined by the selected solcast file
timeInterval = 60/PTM;
nrMinDays = 3; % Number of days over which minima/maxima are averaged
nrMaxDays = 3;
nrAvDays = 3;
timeRes = 15; % time resolution for a profile with user defined time resolution (in minutes, integer)
% Parameters of the PV panels
PVbaseEfficiency = 0.17;
optimalTemperature = 23; 
temperatureCoefficient = -0.05;
panelTemperatureFactor = 1.2;
installedPeakPower = 4; %! in kWp
nrMonths = 24;
goodDayThreshold = 10;
relRiskThreshold = 8; 
midRiskThreshold = 4; 
highRiskThreshold = 2;

%% Save Variables
save SolarDataInputs.mat panelTiltAngle DHI DNI GHI zenith azimuth temperature PTM timeInterval ...
    nrMinDays nrMaxDays nrAvDays PVbaseEfficiency optimalTemperature temperatureCoefficient panelTemperatureFactor ...
    installedPeakPower nrMonths goodDayThreshold relRiskThreshold midRiskThreshold highRiskThreshold timeRes