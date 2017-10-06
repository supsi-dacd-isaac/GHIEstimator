function [GHI_est,GHI_est_i,thetas,err_function,y_hat] = get_ghi(data,identification_pars,optimization_pars,location)
% Estimate GHI from PV power signal and ambient temperature profile.
% Inputs: 
%          data: structure contatining the following fields: {'P','t','T',
%               'GHImcr','GHImcr_clear','GHI3e','GHI'}. All fields but
%               'P','t' and 'T' can be empty (empty array []).
%                P is a T*N matrix, where T is the total number of
%                observations and N is the number of PV plants, containing
%                the power observations.
%                t is a datenum vecto of times relative to the observations
%                T is a vector containing ambient temperature observations
%                GHImcr,GHI3e,GHI are optionally used for comparison, the
%                first two are GHI signals from satellite based services
%                and the last one the GHI signal from a local pyranometer.
%                GHImcr_clear is a satellite-corrected clear sky model,
%                this is also optional. If not provided, a simpler clear
%                sky model is used. 
% 
%         identification_pars: structure containing the following fields: 
%                       {'clearSkyFromSat','n_clear_window','equally_spaced'}. 
%                       clearSkyFromSat is a boolean, if true and if 
%                       GHImcr_clear has been provided in the data
%                       structure, the satellite-corrected clear sky model
%                       is used for identification. n_clear_window specify
%                       how many consecutive samples must be identified as
%                       clear sky condition to consider a period as clear
%                       sky period. equally_spaced is boolean, if true the
%                       proxies are sampled to have a uniform orientation
%                       on a unit sphere, otherwise their orientation is
%                       sampled from a rectangular grid of sun azimuth and
%                       elevation.
% 
%         optimization_pars: sructure containing the following fields:
%         {'expWeight','do_iter_plots','max_iter','lambda',
%          'k_outlier','use_ftrust','n_init','use_single_pv_plant',
%          'correction'}.
%           do_iter_plots is a boolean triggering some plots during 
%           estimation process. 
%           max_iter is a scalar indicating the maximum number of iteration
%           of the optimization algorithm. 
%           lambda is the initial step size of the gradient descent
%           algorithm. 
%           k_outlier set the parameter for the interquantile outlier
%           detection. 
%           use_ftrust is a boolean, if true the trust function is used for
%           estimating GHI.
%           n_init: number of initial search points for the GHI grid search
%           use_single_pv_plant: if true the GHI is estimated N more times
%           using only one PV power plant signal at the same time.
%           correnction: a string in {'T','TL','TLV'}. 'T' applies a
%           temperature correction, 'TL' applies an additional correction
%           considering th inverter efficiency, 'TLV' consider also the
%           influence of low incidence angles.
%          
%         location: a structure containing the following fields: {'UTC',
%         'latitude','longitude','altitude','n_thetas'}. UTC is a scalar 
%         containing the UTC shift in hours, latitude and longitude are 
%         floats containing GPS coordinates and altitude is float 
%         containing the altitude in meters.
% Outputs: 
%          GHI_est: GHI estiamted with all the available signals
% 
%          GHI_est_i: T*N matrix of GHI estimated using only one PV power 
%                     plant signal at the same time.
% 
%          thetas: N_p*N, where N_p is the number of proxies. Matrix of 
%                  identified coefficients relative to the proxies. 
% 
%          err_function: cell array of error functions for the PV signals. 
%                        This is a function of sun azimuth and elevation.
% 
%          y_hat: T*N matrix containing the estimated generated power of
%          the PV plants, given the estimated GHI and estimated PV plant
%          nominal power and orientations.
%          
%% Check presence of all inputs
data_fields = {'P','t','T','GHImcr','GHImcr_clear','GHI3e','GHI'};
identification_fields = {'clearSkyFromSat','n_clear_window','equally_spaced'};
optimization_fields = {'expWeight','do_iter_plots','max_iter','lambda','k_outlier','use_ftrust','n_init','use_single_pv_plant','correction'};
location_fields = {'latitude','longitude','altitude','UTC'};
try
    P = data.P;
    t = data.t;
    T = data.T;
    GHI_sat_clear = data.GHImcr_clear; % can be empty. Used for identification. If not provided, a clear sky model is used instead.
    GHI_sat = data.GHImcr; % can be empty, used for valdiation only
    GHI_sat_2 = data.GHI3e; % can be empty, used for validation only
    GHI = data.GHI; % can be empty, in this case do not plot performances
catch
    for i=1:length(data_fields)
        missing_dat(i) = ~isfield(data,data_fields{i});
    end
    fprintf('Missing data fields:\n')
    sprintf('%s\n',data_fields{missing_dat})
    error('Missing fields in data struct')
end


for i=1:length(identification_fields)
    missing_id(i) = ~isfield(identification_pars,identification_fields{i});
end

if any(missing_id)
    fprintf('Missing identification fields:\n')
    sprintf('%s\n',identification_fields{missing_id})
    error('Missing fields in identification struct')
end

try
     expWeight = optimization_pars.expWeight;
     do_plots = optimization_pars.do_iter_plots;
     lambda = optimization_pars.lambda;
     k_outlier = optimization_pars.k_outlier;
     use_ftrust = optimization_pars.use_ftrust;
     n_init = optimization_pars.n_init;
     use_single_pv_plant = optimization_pars.use_single_pv_plant;
     correction = optimization_pars.correction;
     max_iter = optimization_pars.max_iter;
catch
    for i=1:length(optimization_fields)
        missing_opt(i) = ~isfield(optimization_pars,optimization_fields{i});
    end
    fprintf('Missing optimization fields:\n')
    sprintf('%s\n',optimization_fields{missing_opt})
    error('Missing fields in optimization struct')
end

 
for i=1:length(location_fields)
    missing_loc(i) = ~isfield(location,location_fields{i});
end
if any(missing_loc)
    fprintf('Missing optimization fields:\n')
    sprintf('%s\n',location_fields{missing_loc})
    error('Missing fields in location struct')
end


%% Identify nominal powers and orientations
fprintf('\n-------- Starting PV fields identification --------')
[thetas,err_function] = get_thetas(P,t,location,T,identification_pars,GHI_sat_clear);
thetas = thetas.thetas_biweight; 

%% Estimate GHI
fprintf('\n-------- Starting GHI estimation --------')
opt_pars.expWeight = expWeight;
opt_pars.do_plots = do_plots;
opt_pars.lambda = lambda;
opt_pars.max_iter = max_iter;
opt_pars.correction = correction;
opt_pars.equally_spaced = identification_pars.equally_spaced;

% Use normalized power signals and normalized thetas
P_norm = P./sum(thetas);
thetas_norm =thetas./sum(thetas);

% Estimate GHI recursively
if use_ftrust
    [GHI_est,y_hat] = ghi_estimator(P_norm,GHI,t,location,thetas_norm,T,n_init,opt_pars,k_outlier,err_function);
else
    [GHI_est,y_hat] = ghi_estimator(P_norm,GHI,t,location,thetas_norm,T,n_init,opt_pars,k_outlier,[]);
end

%de-normalize estimated PV signals
y_hat = y_hat.*sum(thetas); 

if use_single_pv_plant
    for i=1:size(P,2)
        P_norm_i = P(:,i)/sum(thetas(:,i));
        thetas_norm_i =thetas(:,i)/sum(thetas(:,i));
        [GHI_est_i(:,i)] = ghi_estimator(P_norm_i,GHI,t,location,thetas_norm_i,T,n_init,opt_pars,k_outlier); 
    end
else
    GHI_est_i = [];
end
    
%% Do plots
if ~isempty(GHI)
    saveplot = false;
    title_str = [];
    do_some_plots(GHI,GHI_est,t,GHI_sat,GHI_sat_2,saveplot,title_str)
end


