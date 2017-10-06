function [error_functions] = map_pv_errors(P,P_est,t,location,threshold,quant)
% Retrieve errors and relative errors
% 
% Inputs:P: T*N matrix of PV plants power observations, where T is the
%             total number of observations and N is the number of PV power
%             plants.
%       
%        P_est: T*N matrix of PV plants simulated power observations,
%        generated with clear GHI signal.
% 
%        t: datenum time vector associated with the P observations.
% 
%        location: a structure containing the following fields: {'UTC',
%         'latitude','longitude','altitude','n_thetas'}. UTC is a scalar 
%         containing the UTC shift in hours, latitude and longitude are 
%         floats containing GPS coordinates and altitude is float 
%         containing the altitude in meters. 
% 
%        threshold: threshold for the error_function.
%       
%        quant: quantile of the error distribution considered as clear
%        observation (eq. 24 in the paper)
% 
% Outputs: err_functions: cell array of error functions for the PV signals. 
%                         This is a function of sun azimuth and elevation.


if nargin<6
    quant=0.01;
end
if nargin<5
    threshold = 0;
end

err_PV = P-P_est;
rel_err_PV = abs(P-P_est)./P;

UTC = location.UTC;
L_lat = location.latitude;
L_long = location.longitude;
L_alt = location.altitude;


%% Compute day of the year. Do not use days365, is deadly slow
[yy,~] = datevec(t);
dayOfYear = floor(t-datenum(yy,1,1))+1;

%% Get Sun position

Time = pvl_maketimestruct(t, ones(size(t))*UTC);            % generate a structured Time
Location = pvl_makelocationstruct(L_lat,L_long,L_alt);      % generate a structured Location
[SunAz, SunEl, AppSunEl, ~] = pvl_ephemeris(Time,Location); % compute azimuth, elevation (90�-zenit) , apparent elevation and solar time. NB:can take additional arguments
Zenit=90-SunEl;
HExtra = pvl_extraradiation(dayOfYear);                     % extraterrestrial irradiance
SunZen = 90 - AppSunEl;                                     % real zenith
AM = pvl_relativeairmass(SunZen);                           % air mass
AM(isnan(AM)) = 20;

%% Get Error matrices
res=5; %deg integer
ix1=round(SunAz/res)+1;
ix2=round(SunEl/res)+1;
filter = SunEl>0;
[Xobs,Yobs] = meshgrid(0:res:360,0:res:90);
for i=1:size(err_PV,2)
    filter = filter & isfinite(rel_err_PV);
    RE_i(:,:,i) = accumarray([ix1(filter),ix2(filter)],rel_err_PV((filter),i),[360/res+1 90/res+1],@(x) quantile(x,quant),nan);
    num_i(:,:,i) = accumarray([ix1(filter),ix2(filter)],rel_err_PV((filter),i),[360/res+1 90/res+1],@(x) length(x),0);    % number of observations per bin
end


%% Plot PV errors as function of SunZen and SunAz
[X,Y] = meshgrid(0:res:360,0:res:90);
quant_num_threshold = 0.3;
for i=1:size(P_est,2)
    RE_dummy = RE_i(:,:,i)';
    num_dummy = num_i(:,:,i);
    only_these = isfinite(RE_dummy);
    RE_dummy(num_i(:,:,i)' < quantile(num_dummy(num_dummy>0),quant_num_threshold))=0;
    RF{i} = fitrgp([X(isfinite(RE_dummy)),Y(isfinite(RE_dummy))], RE_dummy(isfinite(RE_dummy)),'sigma',0.01);
    error_functions{i} = @(az,el) trust_threshold(RF{i},az,el,threshold);
end

res=1;
[X,Y] = meshgrid(0:res:360,0:res:90);
for i=1:size(P_est,2)
    figure
    num_dummy = num_i(:,:,i);
    only_these = num_i(:,:,i)'>quantile(num_dummy(num_dummy>0),quant_num_threshold);
    RE_dummy = RE_i(:,:,i)';
    num_dummy = num_i(:,:,i)';
    % trusted = reshape(predict(RF{i},[X(:),Y(:)]),size(X,1),[]);
    trusted = reshape(error_functions{i}(X(:),Y(:)),size(X,1),[]);
    trusted(abs(trusted)>10) = 0;
    nanMask=imresize(isnan(RE_dummy),size(X),'bilinear')>0.5;
    trusted(nanMask) = NaN;
    contourf(X,Y,trusted,10,'LineStyle','none');hold on;
    xlabel('sunAz (deg)')
    ylabel('sunEl (deg)')
    xlim([50 310])
    ylim([0 70])
    h=colorbar;
    cm=parula;
    cm=flipud(cm);
    colormap(cm)
    h=xlabel(h,'thresholded relative error');
    set(h,'fontsize',7)
end
end

function trusted = trust_threshold(RF,az,el,threshold)
if isa(RF,'RegressionGP')
    trusted = predict(RF,[az,el]);
else
    trusted =  RF(az,el);
end
trusted(trusted<threshold) = 0;
end

