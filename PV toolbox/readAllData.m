function [DHI, DNI, GHI, zenith, azimuth,temperature, consumption] = readAllData(...
    solcastFile, consumptionFile, DHIRange, DNIRange, GHIRange, ...
    zenithRange, azimuthRange, temperatureRange, consumptionRange)
% READALLDATA Reads the variables from an excel provided by Solcast and the
% consumption data
%   Assumes all provided arrays have the same length
%
%   Parameters:
%       - file String with filename to read data
%       - nbParams: Integer, number of parameters to be read
%       - DHICols: 1D array with string indicating what row range in what
%                   column to read from. Positions in array match with
%                   previous parameters
%       - DNICols: 1D array with string indicating what row range in what
%                   column to read from. Positions in array match with
%                   previous parameters
%       - GHICols: 1D array with string indicating what row range in what
%                   column to read from. Positions in array match with
%                   previous parameters
%       - zenithCols: 1D array with string indicating what row range in 
%                   what column to read from. Positions in array match with
%                   previous parameters
%       - azimuthCols: 1D array with string indicating what row range in 
%                   what column to read from. Positions in array match with
%                   previous parameters
%       - temperatureCols: 1D array with string indicating what row range in 
%                   what column to read from. Positions in array match with
%                   previous parameters
%   Returns:
%       - DHI: Diffuse horizontal irradiance. 1D float array of size N x 1
%       - DNI: Direct normal irradiance. 1D float array of size N x 1
%       - GHI: Global horrizontal irradiance. 1D float array of size N x 1
%       - zenith: 1D float array of size N x 1
%       - azimuth: 1D float array of size N x 1
%       - temperature: 1D float array or size N x 1

DHI = []; DNI = []; GHI = []; zenith = []; azimuth = []; temperature = [];
consumption = [];

% Diffuse horizontal irradiance
DHI = [DHI; xlsread(solcastFile,DHIRange)];

% Direct normal irradiance
DNI = [DNI; xlsread(solcastFile,DNIRange)];

% Global horrizontal irradiance
GHI = [GHI; xlsread(solcastFile,GHIRange)];

% Zenith
zenith = [zenith; xlsread(solcastFile,zenithRange)];

% Azimuth
azimuth = [azimuth; xlsread(solcastFile, azimuthRange)];

% Temperature
temperature = [temperature; xlsread(solcastFile,temperatureRange)];

% Consumption
consumption = [consumption; xlsread(consumptionFile,consumptionRange)];

end