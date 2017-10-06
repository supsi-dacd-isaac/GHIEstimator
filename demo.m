% This script runs the GHI estimation algorithm presented in [1] on a 
% partial version of the paper dataset. Only freely available satellite
% data is used. 
% [1] "L.Nespoli, V.Medici, An Unsupervised Method for 
% Estimating the Global Horizontal Irradiance from Photovoltaic Power 
% Measurements, Solar Energy 2017"

clear
clc 
close all

% Load the dataset

% add paths
addpath('src/PV_lib_132'); % Change this with the path of installation of the PV lib for matlab 
addpath(genpath('src'))

% load data
load('data/dataset.mat')

% Define parameters for the GHI estimation
location.latitude = 47.50685;
location.longitude = 7.519069;
location.altitude = 323;
location.UTC = 0;

opt_pars.expWeight = false;
opt_pars.do_iter_plots = false;
opt_pars.max_iter = 30;
opt_pars.lambda = 1;
opt_pars.k_outlier = 1.5;
opt_pars.use_ftrust = true;
opt_pars.n_init = 30;
opt_pars.use_single_pv_plant = false;
opt_pars.correction = 'TLV'; % correct for temperature, low invertert powers and for low irradiation angles

id_pars.clearSkyFromSat = true;
id_pars.n_clear_window = 3;
id_pars.equally_spaced = true; 

% estimate PV models and GHI
[GHI_est,GHI_est_i,thetas,ftrust,y_hat] = get_ghi(data,id_pars,opt_pars,location);
