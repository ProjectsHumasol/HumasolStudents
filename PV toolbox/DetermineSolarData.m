% This file computes information relative to the solar irradiation and the
% PV production for given timeseries.
% The functions used are part of the SolarAnalyzer class, but can also be
% used alone (as they are "static") 
% All the functions mentioned here were designed for the dimensioning of an
% installation using historical data from Solcast: https://solcast.com/historical-and-tmy/?gclid=CjwKCAjwzaSLBhBJEiwAJSRokmsCZHSpSvPPHijRgICzxID4ykgWLe-XVB75X_KaAS6VDWoGck8XKhoC-mQQAvD_BwE
% As a student you can get free credits allowing to request the full
% timeseries for a given location. This should then be used to allow a more
% accurate dimensioning of the installation than with typical "Equivalent
% PV hours"

% This file is meant as an example showing the information that can be
% retrieved from the timeseries using the SolarAnalyzer class. The input
% data should be modified according to your project location. Extensive
% documentation for all functions is present in the SolarAnalyzer class
% More information can be found in the ReadMe file on GitHub
close all
clear all
clc
% First create a structure allowing to store all the information computed
SolarData = struct();

% Load user defined parameters for the computation of PV information
% These should be stored in a mat file and the name below should match the
% name of your file
load("SolarDataInputs.mat") % An example file to create this is also given: "InputsSolarData"

%% Compute a number of auxiliary parameters for the orientation and location
% of the solar panels

% Sky view factor (Parameter determining quantity of irradiation from
% reflection absorbed by a solar panel)
SolarData.SVF = SolarAnalyzer.computeSkyViewFactor(panelTiltAngle);

% Optimal Azimuth (angle between the panels and the line due north), this
% is basically the average between sunrise and sunset
SolarData.OptimalAzimuth = SolarAnalyzer.computeOptimalPanelAzimuth(azimuth, GHI);

%% Compute the solar irradiation for all considered timesteps
% Total solar irradiation computed from the solar data for the
% defined PV configuration
SolarData.totalIrradiation = SolarAnalyzer.irradiation(GHI, DHI, DNI, azimuth,...
                            zenith, panelTiltAngle);
                        
                                    
% Compute average PV efficiency for a given type of solar panel, taking
% into account the efficiency variations due to temperature
SolarData.PVefficiency = SolarAnalyzer.pvPanelEfficiency(PVbaseEfficiency, optimalTemperature,...
                temperatureCoefficient, temperature, panelTemperatureFactor, GHI);
            
% Use the irradiation and panel efficiency to compute the expected PV profile           
[SolarData.sunProfile, SolarData.pvProfile, SolarData.pvProfileRes,SolarData.totalYield, SolarData.averageMonthlyYield,...
                    SolarData.pvHours, SolarData.pvDays, SolarData.good, SolarData.relRisk,...
                    SolarData.midRisk, SolarData.highRisk]...
                        = SolarAnalyzer.computeSolarProfiles(SolarData.totalIrradiation, timeInterval,...
                        installedPeakPower,nrMonths, goodDayThreshold, relRiskThreshold, midRiskThreshold,...
                            highRiskThreshold,timeRes);



% Compute averages, minima and maxima over the whole timeseries 
% This returns the minimal/maximal/average solar production (for a whole day and
% a profile for the minimal day)
[SolarData.minProfile,SolarData.minDay] = SolarAnalyzer.minDailyPV(...
                                        SolarData.pvProfile, 60, nrMinDays);
                                    
[SolarData.avProfile, SolarData.avDay] = SolarAnalyzer.avgDailyPV(...
                                        SolarData.pvProfile, 60);
                                    
[SolarData.maxProfile, SolarData.maxDay] = SolarAnalyzer.maxDailyPV(...
                                        SolarData.pvProfile, 60, nrMaxDays);

% Plot profiles for the minimal, average and maximal irradiation over the
% considered period
SolarAnalyzer.plot(SolarData.minProfile,SolarData.avProfile,SolarData.maxProfile,"2015-2016",60)
% Plot a histogram showing the distribution of PV production over the year
SolarAnalyzer.computePVHistogram(SolarData.pvProfile,installedPeakPower)

if readConsumption 
    % This will only be executed if you change readConsumption
    % to true in the InputsSolarData file.
    % You need an excel file with the consumption data.
    
    % Calculates the battery profile, the time indices giving black-outs
    % and the total black-out time for a specific peakPower and battery
    % capacity.
    [battery, black_outs, black_out_time] = SolarAnalyzer.simulateBattery(...
    	batteryCapacity, batteryEfficiency, minimal_battery_charge, SolarData.pvProfileRes, ...
        consumption, PTM);

    % Calculates the pareto boundary (every combination of battery
    % capacity and peak power that is not dominated by cost, amount of
    % black-outs and black-out time (the last two criteria are very
    % similar, but not exactly the same).
    pareto = SolarAnalyzer.bestCombos(SolarData.pvProfileRes, installedPeakPower, ...
        minPeakPower, maxPeakPower, peakPowerStep, minBatteryCap, ...
        maxBatteryCap, batteryCapStep, batteryEfficiency, minimal_battery_charge, ...
        consumption, PTM, cost);
    
    % Plot battery in certain period
    SolarAnalyzer.plotBattery(battery, '1-mar-2015', '31-mar-2015', '1-jan-2015', PTM)
    % Plot pareto combo (for black-out time, not for total amount of
    % black-outs)
    SolarAnalyzer.plotPareto(pareto)
   
    if optimizeOrientation
        % This calculates the pareto boundary while optimizing panel
        % orientation (how many panels south, east and west).
        paretoVariableRotationAngle = SolarAnalyzer.bestCombosVariableRotationAngle(...
            minPeakPower, maxPeakPower, peakPowerStep, minBatteryCap, ...
            maxBatteryCap, batteryCapStep, batteryEfficiency, minimal_battery_charge, ...
            consumption, PTM, cost, GHI, DHI, DNI, azimuth, zenith, ...
            panelTiltAngle, 0.1, timeInterval, installedPeakPower,nrMonths, ...
            goodDayThreshold, relRiskThreshold, midRiskThreshold,...
            highRiskThreshold,timeRes);
    
        SolarAnalyzer.plotPareto(paretoVariableRotationAngle)
    end
end

save("SolarDataExample.mat","SolarData")