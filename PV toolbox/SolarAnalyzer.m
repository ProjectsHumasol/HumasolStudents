classdef SolarAnalyzer
    % SOLARANALYZER Class containing static methods to analyze and prepare
    %       solar data
    %   This class contains static methods that can be useful for computing
    %   the solar irradiance, as well as interesting statistics about it.
    %   It can also be used to further prepare the data for a simulation of
    %   an electrical PV system
    
    methods
        function obj = SolarAnalyzer()
            % SOLARANALYZER Construct an instance of this class
            
        end
    end
    
    methods (Static)
        function svf = computeSkyViewFactor(panelTiltAngle)
            % COMPUTESKYVIEWFACTOR Computes the sky view factor of the
            %       panels in function of the panel tilt angle
            %   Parameters:
            %       - panleTiltAngle: tilt angle of the solar panels.
            %                   Float, in degrees
            %   Returns:
            %       - svf: Float. Sky view factor
            
            svf = (1 + cosd(panelTiltAngle))/2;
        end
        
        function optimalAzimuth = computeOptimalPanelAzimuth(azimuthData, GHIData)
            % OPTIMALAZIMUTH Computes the optimal azimuth angle for the
            %       solar panels (angle between the panels and the line due
            %       north).
            %   Computed the average azimuth angle between the extremes form
            %   the provided data array. Only the angles when the sun is up
            %   are considered. Panels placed with this angle will capture
            %   the most amount of sunlight throughout the day
            %   Parameters:
            %       - azimuthData: 1D float array with solar azimuth
            %                   data
            %       - GHIData: 1D float array with the global horizontal
            %                   irradiation matching the time points in the
            %                   azimuthData array. Must be of the same
            %                   length as that array
            %   Returns:
            %       - optimalAzimuth: float indicating the optimal azimuth
            %                   angle
            
            % Check that input data is of correct size
            assert(min(size(azimuthData)) == 1);
            assert(SolarAnalyzer.sameSize(azimuthData, GHIData));
            
            % Prealocate memory for relevant data
            relevantAzimuth = zeros(nnz(GHIData), 1);
            
            % Loop over all time points
            k = 1;
            for i = 1:length(azimuthData)
                % Only look at azimuth when sun is up
                if GHIData(i) ~= 0
                    relevantAzimuth(k) = azimuthData(i); 
                    k = k + 1;
                end
            end
            
            % Compute optimal angle
            M = max(relevantAzimuth);
            m = min(relevantAzimuth);
            optimalAzimuth = (M + m)/2;
        end
        
        function totalIrradiation = irradiation(GHI, DHI, DNI, azimuth,...
                            zenith, panelTiltAngle, alpha, panelRotationAngle)
            % IRRADIATION Computes total solar irradiation on the solar
            %       panels
            %   
            %   Parameters:
            %       - GHI: 1D float array. Global horizonal irradiance
            %       - DHI: 1D float array. Diffuse horizontal irradiance.
            %                   Must have the same length as the previous
            %                   one
            %       - DNI: 1D float array. Direct normal irradiance. Must
            %                   have the same length as the previous ones
            %       - azimuth: 1D float array. The angle between a line
            %                   pointing due north to the sun's current 
            %                   position in the sky. 
            %                   Negative to the East. Positive to the West.
            %                   0 at due North. Must be in degrees, and
            %                   must have the same length as the previous
            %                   ones
            %       - zenith: 1D float array. The angle between a line 
            %                   perpendicular to the earth's surface and the
            %                   sun (90 deg = sunrise and sunset; 0 deg = sun 
            %                   directly overhead). Must be in degrees and have
            %                   the same length as previous arrays.
            %       - panelTiltAngle (optional): Float value, in degrees.
            %                   Defaults to 6 degrees.
            %                   Solar panel tilt angle, w.r.t. the 
            %                   horizontal
            %       - alpha (optional): Float value, between 0 and 1. 
            %                   Defaults to 0.10. The albedo or 
            %                   Reflection coefficient, depends on
            %               	soil conditions at the solar panel location
            %       - panelRotationAngle (optional): Float value, in degrees.
            %                   Defaults to 0 degrees. Solar panel angle 
            %                   (rotation about vertical).
            %                   Southern hemisphere, use default value: solar 
            %                   panels facing north (0 deg).
            %   Returns:
            %       - totalIrradiation: 1D float array of the same length
            %                   as the input arrays.
            %                   Total solar irradiation absorbed by the
            %                   solar panels.
            
            % Provide default values if necessary
            if (nargin < 6)
                panelTiltAngle = 6; 
                panelRotationAngle = 0; 
                alpha = 0.10; 
            elseif (nargin < 7)
                alpha = 0.10;
                panelRotationAngle = 0; 
            elseif (nargin < 8)
                panelRotationAngle = 0; 
            end
            
            % Check input data
            assert(min(size(GHI)) == 1);
            assert(SolarAnalyzer.sameSize(GHI, DHI, DNI, azimuth, zenith))
            assert(alpha >= 0 && alpha <= 1);
            azimuth = mod(azimuth, 360);
            zenith = mod(zenith, 360);
            panelTiltAngle = mod(panelTiltAngle, 360);
            panelRotationAngle = mod(panelRotationAngle, 360);
            
            
            % Compute the irradiation using diffuse, direct, global
            % irradiation and the relative angles of panels and sun

            SVF = SolarAnalyzer.computeSkyViewFactor(panelTiltAngle);

            sunAltitude = max(90 - zenith, 0);     % max ensures no negative values

            % Initialize empty vectors of right size, for efficiency
            dirIrradiation = zeros(length(GHI),1);
            groundIrradiation = zeros(length(GHI),1);
            difIrradiation = zeros(length(GHI),1);
            totalIrradiation = zeros(length(GHI),1);
            
            % Compute irradiation
            for i=1:length(GHI)
               if sunAltitude(i) > 0       % && AM-90 <= Azimuth(i)  && Azimuth(i) <= AM+90
                   
                    dirIrradiation(i) = DNI(i) * (sind(panelTiltAngle) * ...
                        cosd(sunAltitude(i)) * cosd(panelRotationAngle - azimuth(i)) + ...
                        cosd(panelTiltAngle) * sind(sunAltitude(i)));
                    
                    groundIrradiation(i) = GHI(i) * alpha * (1 - SVF);
                    
                    difIrradiation(i) = DHI(i) * SVF;
                    
                    totalIrradiation(i) = dirIrradiation(i) +...
                                groundIrradiation(i) + difIrradiation(i);
               else
                    totalIrradiation(i) = 0;
               end
            end
        end
        
        function [minInstant, minDay] = minDailyPV(...
                                        PVprofile, spH, nrMinDays)
            % MINDAILYPV Computes the minimal PV production per
            %       time point and hour of the day.
            %   Uses a certain number of days, nrMinDays, of minimal solar
            %   production to average the minimal production per time
            %   sample and per hour.
            %   Parameters:
            %       - PVprofile: 1D float array.
            %                   Total solar power produced by the
            %                   solar panels.
            %       - spH: integer. number of samples per hour (60 for
            %       solar data computed through solarProfile)
            %       - nrMinDays: integer. How many days with minimal
            %                   irradiation to use to average the value of
            %                   the irradiation.
            %   Returns:
            %       - minInstant: 1D float array, size (24 * spH) x 1, 
            %                   of (average) minimal irradiation values per
            %                   time point of the day. Spacing of
            %                   measurments is the same as the provided
            %                   data array.
            %       - minDay: float. Minimal irradiation for a whole day
            %       (average over the selected number of days)
            
            samplesPerDay = 24 * spH;
            lengthData = length(PVprofile);
            nrDays = floor(lengthData/samplesPerDay);
            
            days = zeros(nrDays, 1);
            
            % Compute total solar irradiation per day
            for i = 1:nrDays
                days(i) = sum(PVprofile( ((i-1) * samplesPerDay + 1) : (i * samplesPerDay)) );
            end
            
            % Find nrMinDays days with minimal irradiation
            [~,lowDays] = mink(days,nrMinDays);         
            lowDays = lowDays - 1;                      % Indexing in matlab starts from 1

            % Calculate total irradiation per time point for these days
            minInstant = zeros(samplesPerDay,1);
            for k = 1:samplesPerDay
                minInstant(k) = sum(PVprofile(lowDays * samplesPerDay + k));
            end
            
            % Average per time point over selected days
            minInstant = minInstant/length(lowDays);
            % Average per hour
            minDay = sum(minInstant)/spH; 
        end
        
        function [avgInstant, avgDay] = avgDailyPV(...
                                        PVprofile, spH)
            % AVGDAILYIRRADIATION Computes average PV production for every
            %       time point of the day for the given data.
            %   Based on https://web.stanford.edu/group/efmh/jacobson/Articles/I/TiltAngles.pdf
            %   Same result with https://arrow.dit.ie/cgi/viewcontent.cgi?article=1030&context=aaconmusart
            %   Parameters:
            %       - PVprofile: 1D float array.
            %                   Total solar power produced by the
            %                   solar panels.
            %       - spH:  integer. number of samples per hour (60 for
            %                   solar data computed through solarProfile)
            %   Returns:
            %       - avgInstant: 1D float array, size (24 * 60/PTM) x 1,
            %                   of average total production per measurment
            %                   time point.
            %       - avgDay: float. float. Average production for a whole day
            %       (average over the selected number of days)
            
            samplesPerDay = 24 * spH;
            lengthData = size(PVprofile,1);
            nrDays = floor(lengthData/samplesPerDay);
            
            avgInstant = zeros(samplesPerDay,1);
            
            % Compute total irradiation collected for each time point (over
            % all days)
            for k = 1:samplesPerDay
                avgInstant(k) = sum(PVprofile((0:nrDays-1) * samplesPerDay + k));
            end
            
            % Compute average
            avgInstant = avgInstant/nrDays;
            avgDay = sum(avgInstant)/spH;
        end
        
        function [maxInstant, maxDay] = maxDailyPV(...
                                        PVprofile, spH, nrMaxDays)
            % MAXDAILYPV Computes the maximum pv production per time
            %       point of the day.
            %   Uses a certain number of days, nrMaxDays, to average the
            %   maximal solar power collected.
            %   Parameters:
            %       - PVprofile: 1D float array.
            %                   Total solar power produced by the
            %                   solar panels.
            %       - spH:  integer. number of samples per hour (60 for
            %                   solar data computed through solarProfile)
            %       - nrMaxDays: integer. How many days with maximal
            %                   irradiation to use to average the value of
            %                   the irradiation.
            %   Returns:
            %       - maxInstant: 1D float array, size (24 * 60/PTM) x 1, 
            %                   of (average) maximal irradiation values per
            %                   time point of the day. Spacing of
            %                   measurments is the same as the provided
            %                   data array.
            %       - maxDay: float. float. Minimal irradiation for a whole day
            %       (average over the selected number of days)
            
            samplesPerDay = 24 * spH;
            lengthData = length(PVprofile);
            nrDays = floor(lengthData/samplesPerDay);
            
            days = zeros(nrDays, 1);
            
            % Find days with most irradiation
            for i = 1:nrDays
                days(i) = sum(PVprofile( ((i-1) * samplesPerDay + 1):(i * samplesPerDay) ));
            end
            
            [~,highDays] = maxk(days,nrMaxDays);
            highDays = highDays - 1;        % Indexing in matlab starts from 1
            
            maxInstant = zeros(samplesPerDay, 1);
            % Calculate total irradiation per time point for these days
            for k = 1:samplesPerDay
                maxInstant(k) = sum(PVprofile(highDays * samplesPerDay + k));
            end
            
            % Average irradiation
            maxInstant = maxInstant/length(highDays);
            maxDay = sum(maxInstant)/spH;
        end
        
        function plot(min, avg, max, period,PTM, isSubfigure, time)
            % PLOT Plots the minimum, average and maximum irradiation over
            %       the given time span.
            %   This function can be used to plot a single figure, or as
            %   subplots for the same figure.
            %   Parameters:
            %       - time: 1D float array of time points to use on the x
            %                   axis.
            %       - min: 1D float array of minimal irradiation values.
            %                   Must have the same size as the time array.
            %       - avg: 1D float array of average irradiation values.
            %                   Must have the same size as the time array.
            %       - max: 1D float array of maximum irradiation values.
            %                   Must have the same size as the time array.
            %       - period: string. Name of the period to be plotted.
            %                   Used in the title of the plot.
            %       - isSubfigure (optional): boolean. Indicates whether 
            %                   this plot should be part of a subfigure or 
            %                   create its own figure. Defaults to false.
            
            if nargin < 6
                nbS = size(max);
                nbS = nbS(1);
                time = linspace(1,nbS/PTM,nbS);
                isSubfigure = 0;
            elseif nargin <7
                nbS = size(max);
                nbS = nbS(1);
                time = linspace(1:nbS/PTM,nbS);
            end
            
            if isSubfigure
                plot(time, max, 'color', 'red')
                plot(time, min, 'color', 'blue')
                plot(time, avg, 'color', 'yellow')
                hold off
                title('Irradiation ' + period)
                xlabel('Hours')
                ylabel('Total irradiation [kW]')
                legend("max","min","avg")
            else
                figure
                hold on
                plot(time, max, 'color', 'red')
                plot(time, min, 'color', 'blue')
                plot(time, avg, 'color', 'yellow')
                hold off
                title('Irradiation ' + period)
                xlabel('Hours')
                ylabel('Total irradiation [kW]')
                legend("max","min","avg")
                hold off
            end
        end
        
        function computePVHistogram(pvProfile,peakPower)
            nbS = size(pvProfile);
            nbS = nbS(1);
            houraverages = zeros(nbS/60,1);

            % Compute average irradiation per hour
            for i = 1:nbS/60
                houraverage = 0;
                for j = 1:60
                    houraverage = houraverage + pvProfile(60*i+j-60);
                end
                houraverage = houraverage/(60); % in kWh
                houraverages(i) = houraverage;
            end
            
            % Using these hourly averages compute total PV per day
            days = zeros(nbS/(24*60),1);
            for i = 1:nbS/(60*24)
                day = 0;
                for j=1:24
                    day = day+houraverages(i*24-24+j);
                end
                days(i) = day;
            end
            days = days; %kWh
            % Plot PV production
            figure
            hold on
            plot(days,'color','red')
            title('PV production selected years [kWh/day]')
            xlabel('days')
            ylabel('Total PV production')
            hold off
            % How many days in a year with what production?
            figure
            hold on
            histogram(days,60)
            title('Histogram selected years')
            xlabel('Total daily PV production in kWh')
            ylabel('Number of days')
            hold off

            % How many consecutive days of low production?
            MAX_Low = 1;
            Low = 1 ;
            Consecutive_days = zeros(ceil(peakPower)*60,1);
            for l = 1:ceil(peakPower)*60
                for i = 1:(nbS/(60*24)-1)
                    if days(i)<l/10
                        for k= 1:(731-i)
                            if days(i+k)>l/10
                                break
                            else
                                Low = Low +1;
                            end
                        end
                    end
                    if Low > MAX_Low
                        MAX_Low = Low;
                    end
                    Low = 1;
                end
                Consecutive_days(l)= MAX_Low;
                MAX_Low=1;
            end

            figure
            hold on
            plot(linspace(1,ceil(peakPower)*6,ceil(peakPower)*60),Consecutive_days)
            title('Consecutive days less then X')
            xlabel('Threshold (kWh/day)')
            ylabel('Number of days')
            hold off
            
        end
        
        function effVar = computeEfficiencyDifference(temperatureData,...
                                                    GHIData, Topt, Tcoeff)
            % COMPUTEEFFICIENCYDIFFERENCE Computes the efficiency
            %       difference of the solar panels w.r.t. ideal conditions
            %   Calculates the efficiency loss due to temperature, uses
            %   weighted average with irradiation data as weight.
            %   More accurate values can be obtained by computing the solar
            %   panel temperature and not the ambient temperature as solar
            %   panels will be directly in the sun and thus receive
            %   additional heat.
            %   https://www.civicsolar.com/article/how-does-heat-affect-solar-panel-efficiencies
            % 
            %   Parameters:
            %       - temperatureData: 1D float array. Temperature of the
            %                   solar panels (not ambient temperature) for
            %                   each time point.
            %       - GHIData: 1D float array, size should match the 
            %                   temperature array. Global horizonal
            %                   irradiance.
            %       - Topt: float. Optimal temperature of the solar panels.
            %                   Units should match those for the
            %                   temperature array.
            %       - Tcoeff: float. Temperature coefficient of the solar
            %                   panels. Units should match the temperature
            %                   array and the GHI array.
            %   Returns:
            %       - effVar: float. Efficiency variation of the solar
            %                   panels, given as a percentage.
            
            efficiencyVar = zeros(length(GHIData), 1);
            for i = 1:length(GHIData)
                % Formula to compute efficiency ifo temperature
                efficiencyVar(i) = GHIData(i) * (temperatureData(i) - Topt) * Tcoeff; 
            end
            effVar = sum(efficiencyVar)/sum(GHIData);
        end
        
        
        function efficiency = pvPanelEfficiency(baseEfficiency, optimalTemperature,...
                temperatureCoefficient, airTemperature, panelTemperatureFactor, GHIData)
            % PVPANELEFFICIENCY Computes the efficiency of the solar panels
            %   The solar panels have a predefined efficiency, but the
            %   actual efficiency will vary depending on the operating
            %   conditions. Uses the method computeEfficiencyDifference to
            %   account for the variation in working conditions.
            %   Parameters:
            %       - baseEfficiency: float. Predefined efficiency of the
            %                   solar panels.
            %       - optimalTemperature: float. Optimal working
            %                   temperature of the solar panels. From
            %                   manufacturer, see model sheets.
            %       - temperatureCoefficient: float. Coefficient indicating
            %                   the relationship between the temperature
            %                   and the efficiency. From manufacturer, see 
            %                   model sheets. Used to compute the
            %                   efficiency variation. Units should match
            %                   the temperature data.
            %       - airTemperature: 1D float array. Air temperature for
            %                   each time point where the efficiency should
            %                   be computed.
            %       - panelTemperatureFactor: float. Factor indicating the
            %                   relationship between the ambient
            %                   temperature and the solar panel
            %                   temperature. Provides a more accurate
            %                   efficiency variation. Units should match
            %                   those of the temperature data.
            %       - GHIData: 1D float array, size should match that of
            %                   the air temperature data. Units should
            %                   match the temperature data and the
            %                   temperature coefficient.
            %   Returns:
            %       - efficiency: float. True solar panel efficiency for
            %                   the provided working conditions.
            
            assert(min(size(airTemperature)) == 1 &&...
                SolarAnalyzer.sameSize(airTemperature,GHIData));
            
            tempData = airTemperature * panelTemperatureFactor;
            
            efficiencyvar = SolarAnalyzer.computeEfficiencyDifference(...
                                tempData,GHIData,optimalTemperature,...
                                temperatureCoefficient);
            
            efficiency = baseEfficiency + efficiencyvar/100; % efficiencyvar is returned as percentage
        end
        
        function [sunProfile, pvProfile, pvProfileRes, totalYield, averageMonthlyYield,...
                    pvHours, pvDays, good, relRisk, midRisk, highRisk]...
                        = computeSolarProfiles(solarData, timeInterval, installedPeakPower,...
                            nrMonths, goodDayThreshold, relRiskThreshold, midRiskThreshold,...
                            highRiskThreshold,timeRes)
            % COMPUTESOLARPROFILES Computes the solar production profile an
            %       the risk level per day
            %   Converts the given solar data array into a more
            %   fine-grained data array, data is now available every
            %   minute. Does this by expanding each given value over the
            %   given timeInterval.
            %   Parameters:
            %       - solarData: 1D float array. Total solar irradiance
            %                   collected by the solar panels.
            %       - timeInterval: integer, in minutes. Indicates time
            %                   between measurments in the solarData array.
            %                   Must be positive.
            %       - installedPeakPower: float, in W. Peak power of the
            %                   installed solar panel array.
            %       - nrMonths: integer. Number of months present in the
            %                   data. Used to compute a monthly average.
            %       - goodDayThreshold: float, in Wh. Energy threshold to
            %                   determining if there was enough production
            %                   during a day. Can also be used to see if
            %                   there is eccess productuin (to power other
            %                   loads).
            %       - relRiskThreshold: float, in Wh. Energy threshold to
            %                   indicate when the production is not enough
            %                   for 24 hours.
            %       - midRiskThreshold: float, in Wh. Energy threshold for
            %                   mid risk.
            %       - highRiskThreshold: float, in Wh. Energy threshold to
            %                   indicate when the production is not enough
            %                   if the previous day was risky.
            %   Returns:
            %       - sunProfile: 1D float array of size 
            %                   (length(solarData) * timeInterval) x 1.
            %                   Piecewise extended solar data array. Same
            %                   units as provided data.
            %       - pvProfile: 1D float array of size 
            %                   (length(solarData) * timeInterval) x 1.
            %                   Available solar power from the panels per 
            %                   minute. In W.
            %       - totalYield: float, in Wh. Total solar energy
            %                   production over the time span of the
            %                   provided data.
            %       - averageMonthlyYield: float, in Wh. Average solar
            %                   energy production per month.
            %       - pvHours: 1D float array of size
            %                   (length(solarData) * timeInterval /60) x 1.
            %                   Hourly solar energy production in Wh.
            %       - pvDays: 1D float array of size
            %                   (length(solarData) * timeInterval /(24*60)) x 1.
            %                   Daily solar energy production in Wh.
            %       - good: integer. Count of days with good solar
            %                   production, i.e. that surpass the 
            %                   goodDayThreshold.
            %       - relRisk: integer. Count of days with relative risk,
            %                   i.e. a production below relRiskThreshold.
            %                   Does not exclude middle and high risk days.
            %       - midRisk: integer. Count of days with middle risk,
            %                   i.e. a production below midRiskThreshold.
            %                   Does not exclude high risk days.
            %       - highRisk: integer. Count of days with high risk, i.e.
            %                   a production below highRiskThreshold.
            
            assert(timeInterval > 0)
            
            inputLength = size(solarData, 1);
            outputLength = inputLength * timeInterval;      % In minutes
            pvProfileRes = zeros(outputLength/timeRes,1);
            % Compute solar profile for every minute
            if timeInterval ~= 1
                sunProfile = zeros(outputLength, 1);

                for i = 2:(inputLength - timeInterval)
                    sunProfile(( (i - 1)*timeInterval:i*timeInterval )) =...
                            solarData(i);
                end
            else
                sunProfile = solarData;
            end
            sunProfile = sunProfile/1000; % in kW
            
            
            % PV yield
            pvProfile = sunProfile * installedPeakPower;
            totalYield = sum(pvProfile)/60; % in kWh, divide by 60 because was in Wmin
            averageMonthlyYield = totalYield/nrMonths; % in Wh per month
            
            % Compute PV profile with chosen time resolution
            for i = 1:outputLength/timeRes
                temp = 0;
                for k = 1:timeRes
                    temp = temp+pvProfile((i-1)*timeRes+k);
                end
                pvProfileRes(i) = temp/timeRes;
            end
                
            
            % Compute average irradiation per hour
            pvHours = zeros(floor(outputLength/60),1);
            for i = 1:outputLength/60
                pvHours(i) = sum(pvProfile( (i-1)*60+1:i*60 )) / 60;
            end
            
            % Total PV per day, using hourly averages
            pvDays = zeros(floor(outputLength/(60*24)),1);
            for i = 1:outputLength/(60*24)
                pvDays(i) = sum(pvHours( (i-1)*24+1:i*24 ));
            end
            
            % Check PV yield and risk
            good        = sum(pvDays > goodDayThreshold);
            relRisk     = sum(pvDays < relRiskThreshold);
            midRisk     = sum(pvDays < midRiskThreshold);
            highRisk    = sum(pvDays < highRiskThreshold);
            
        end
        
        function bool = sameSize(varargin)
            % SAMESIZE Checks wheter the provided 2D arrays have the same
            %       size
            %   True if all arrays have the same maximum dimension and the
            %   same minimum dimension. Is usefull for 1D arrays, since it
            %   matters less if it is a row or column vector. Not for truly
            %   2D arrays, since it is not the same as size(x) == size(y).
            %   Parameters:
            %       - varargin: comma sepparated list of 2D arrays.
            %   Returns:
            %       - bool: boolean. True if all provided arrays have the
            %                   same maximum dimension and the same
            %                   minimum dimension, false otherwise.
            
            bool = 1;
            if ~isempty(varargin)
                x = max(size(varargin{1}));
                y = min(size(varargin{1}));
                for k = 1:length(varargin)
                    if ~(max(size(varargin{k})) == x && min(size(varargin{k})) == y)
                        bool = 0;
                        break;
                    end
                end
            end
        end
    end
end

