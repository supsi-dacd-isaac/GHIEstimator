function [thetas,error_functions] = get_thetas(P,t,location,Ta,opt,GHI_clear)
% Identify nominal power and orientations of the PV fields which generated
% the power signal in the P matrix, by means of coefficient of the proxy
% matrix.
% 
% Inputs: P: a T*N matrix, where T is the total number of observations and 
%            N is the number of PV plants, containing the power 
%            observations.
%
%         t: datenum vector relative to observations in P
% 
%         location: a structure containing the following fields: {'UTC',
%         'latitude','longitude','altitude','n_thetas'}. UTC is a scalar 
%         containing the UTC shift in hours, latitude and longitude are 
%         floats containing GPS coordinates and altitude is float 
%         containing the altitude in meters.
%         
%         Ta: vector of ambient temperature observations
% 
%         opt: structure containing the following fields: 
%             {'clearSkyFromSat','n_clear_window','equally_spaced'}. 
%             clearSkyFromSat is a boolean, if true and if  GHImcr_clear 
%             has been provided in the data structure, the 
%             satellite-corrected clear sky model is used for 
%             identification. n_clear_window specify how many consecutive 
%             samples must be identified as clear sky condition to consider
%             a period as clear sky period. equally_spaced is boolean, 
%             if true the proxies are sampled to have a uniform orientation
%             on a unit sphere, otherwise their orientation is sampled from
%             a rectangular grid of sun azimuth and elevation. 
% 
%         GHI_clear: clear sky GHI signal from satellite. Can be empty. In
%         this case the clear GHI signal is obtained from a simpler clear
%         sky model.
% 
% Output: thetas: N_p*N, where N_p is the number of proxies. Matrix of 
%                 identified coefficients relative to the proxies. 
%   
%         err_functions: cell array of error functions for the PV signals. 
%                        This is a function of sun azimuth and elevation.


% fitgmdist initializes with random points, this way one always gets the
% same result
rng(1); 

% retrieve GHI clear sky signal
if nargin<6 || ~opt.clearSkyFromSat
    GHI_clear = ghi_clear_sky(t,location,location.UTC);
end

%% Get clear observations
UTC = location.UTC;
L_lat = location.latitude;
L_long = location.longitude;
L_alt = location.altitude;
Time = pvl_maketimestruct(t, UTC);
Location = pvl_makelocationstruct(L_lat,L_long,L_alt);

 %compute azimuth, elevation (90deg-zenit) , apparent elevation and solar
 %time. NB:can take additional arguments
[SunAz, SunEl] = pvl_ephemeris(Time,Location);
res=5; % size of square window in deg
clear_idx=false(size(P));

for i=1:size(P,2)
    PVSingle=P(:,i);
    filt=SunEl>0&PVSingle>0;
    PVSingleF=PVSingle(filt);
    SunElF=SunEl(filt);
    SunAzF=SunAz(filt);
    testedEl=0:90;
    testedAz=0:360;
    sel=nan(length(testedEl),length(testedAz),2);
    el=nan(length(testedEl),length(testedAz));
    az=nan(length(testedEl),length(testedAz));
    mask=false(length(testedEl),length(testedAz));
    for p=1:length(testedEl)
        for q=1:length(testedAz)
            el(p,q)=testedEl(p);
            az(p,q)=testedAz(q);
            PVSel=PVSingleF(SunElF>(testedEl(p)-res/2)&SunElF<(testedEl(p)+res/2)&SunAzF>(testedAz(q)-res/2)&SunAzF<(testedAz(q)+res/2));
            if ~isempty(PVSel)
                
                % test for bimodality of the distribution. The power
                % distribution tend to be unimodal at low sun elevation
                % angles and bimodal at higher sun elevation angles.
                dip = HartigansDipTest(PVSel);
                
                %if distribution is not bimodal, drop the el,az combination
                if dip <0.02
                    continue;
                end
                
                mask(p,q)=true;
                
                % If the distribution seems bimodal, fit a gaussian mixture
                % on it. 
                try 
                        gm=fitgmdist(PVSel,2,'Options',statset('MaxIter',1e4));
                        [mu,ix]=max(gm.mu);
                        sigma=gm.Sigma(ix);
                        sel(p,q,:)=[mu-sqrt(sigma), mu+sqrt(sigma)];
                end
            end
        end
    end
    
    
    LL=sel(:,:,1);
    toKeep=isfinite(LL(:));
    
    % applying smoothing function (gaussian filter) to lower and upper
    % bounds 
    TK=scatteredInterpolant(az(:),el(:),double(toKeep),'nearest','nearest');
   
    LLS=scatteredInterpolant(az(toKeep),el(toKeep),LL(toKeep),'linear','nearest');
    LLI=reshape(LLS(az,el),length(testedEl),length(testedAz));
    LLG=imgaussfilt(LLI,3);
    LLGP=LLG;
    LLGP(~mask)=nan;
    
    LU=sel(:,:,2);
    toKeep=isfinite(LU(:));
    LUS=scatteredInterpolant(az(toKeep),el(toKeep),LU(toKeep),'linear','nearest');
    LUI=reshape(LUS(az,el),length(testedEl),length(testedAz));
    LUG=imgaussfilt(LUI,3);
    LUGP=LUG;
    LUGP(~mask)=nan;
    
    figure('position',[100,100,1000,400]);
    subplot(121)
    [~,h]=imagescwithnan(LLGP,parula,[1 1 1],0:res:360,0:res:90,[0 max(LUGP(:))]);
    set(gca,'YDir','normal')
    xlabel('azimuth (deg)')
    ylabel('elevation (deg)')
    xlabel(h,'Power (W)')
    title(sprintf('PV %i %s' ,[i,'minimum Power (\mu - \sigma)']))
    subplot(122)
    [~,h]=imagescwithnan(LUGP,parula,[1 1 1],0:res:360,0:res:90,[0 max(LUGP(:))]);
    set(gca,'YDir','normal')
    xlabel('azimuth (deg)')
    ylabel('elevation (deg)')
    xlabel(h,'Power (W)')
    title(sprintf('PV %i %s' ,[i,'maximum Power (\mu + \sigma)']))

    
    % reinterpolate only based on trusted points
    LLS=scatteredInterpolant(az(mask),el(mask),LLG(mask));
    LUS=scatteredInterpolant(az(mask),el(mask),LU(mask));
    
    % identify clear sky periods by checking if the power value is between
    % mu-sigma and mu+sigma and if it is above 0 and below 97% of the maximum power
    % output (this is done to avoid curtailment periods)
    clear_idx(filt,i) =  arrayfun(@(pv,az,el)...
        pv>LLS(az,el)&pv<LUS(az,el)&TK(az,el),PVSingleF,SunAzF,SunElF)&...
        PVSingleF<0.97*max(PVSingleF)&...
        PVSingleF>0;
    
    figure;
    plot(t,PVSingle);
    hold on
    plot(t(clear_idx(:,i)),PVSingle(clear_idx(:,i)),'*')
    % select only periods in which the sky is clear for at least
    % opt.n_clear_window times in a row
    clear_idx(:,i)=find_long_sequences(clear_idx(:,i),opt.n_clear_window);
    plot(t(clear_idx(:,i)),PVSingle(clear_idx(:,i)),'o')
    datetick
    legend({sprintf('PV %i',i),'Clear period','Clear periods with at least n\_clear\_window obs.'})
    ylabel('Power [W]')
    drawnow;
end

%% Get proxies
equally_spaced = opt.equally_spaced;
proxies = compute_proxies(GHI_clear,t,Ta,location,equally_spaced,'TLV');

h = figure;
s = figure;

opt_pars.lambda = 1e-3;
opt_pars.rescale = false;

%% Identify PV fields for each PV power signal and get 
for i=1:size(P,2)
    PV_q(i) = quantile(P(clear_idx(:,i),i),0.9);
    X = proxies(clear_idx(:,i),:);
    target = P(clear_idx(:,i),i);
    filter = target>(0.1*PV_q(i));
    X = X(filter,:);
    target = target(filter);
    
    if opt_pars.rescale
        opt_pars.epsilon = 0.1;
    else
        opt_pars.epsilon = PV_q(i)/20;
    end
    
    t_filt = t(clear_idx(:,i));
    t_filt = t_filt(filter);
    
    % Perform two consecutive robust regressions (the first one uses huber
    % loss and it's used as starting point for the biweight robust
    % regression). This identification procedure has been modified from the
    % one performed in the paper, which used recursive reweighted least
    % squares and was performed with a modified Matlab propietary script.
    % The following implementation uses fmincon solver. 
    
    opt_pars.loss_function = 'huber';
    [thetas_huber(:,i)] = robust_proxy_fit(target,X,opt_pars);
    opt_pars.loss_function = 'biweight';
    [thetas_biweight(:,i)] = robust_proxy_fit(target,X,opt_pars,thetas_huber(:,i));
    
    figure(h)
    subplot(2,2,rem(i-1,4)+1)
    plot(target);hold on;
    plot(X*[thetas_biweight(:,i)],'--');
    legend({'Observations','Simulated'});
    title(sprintf('PV %i',i))
    xlabel('Clear observation number')
    ylabel('Power [W]')
    figure(s)
    subplot(2,2,rem(i-1,4)+1)
    scatter(target,X*thetas_biweight(:,i),'.');
    hold on;
    plot([min(P),max(P)],[min(P),max(P)],'--r','linewidth',2);
    xlabel('Observations [W]')
    ylabel('Simulated [W]')
    title(sprintf('PV %i',i))
    
    drawnow;
    thetas.thetas_biweight = thetas_biweight;
    
    % Build the error function searching for systematic shadows - the error
    % function will be used to build the trust function in ghi_estimator.m
    % function
    threshold = 0.25;
    clear_quantile = 0.01;
    error_functions_i = map_pv_errors(P(:,i),proxies*thetas_biweight(:,i),t,location,threshold,clear_quantile);
    error_functions{i} = error_functions_i{1};
    
end
end