function proxies = compute_proxies(GHI,t,Ta,location,equallySpaced,correction)
% Compute the proxy signals. Proxies are estimation of the power produced
% by a PV panel under different orientations and a given GHI. 
% Inputs: GHI: vector of observations of GHI for which the proxy are computed
%         
%         t: datenum vector relative to the GHI observations
%         
%         Ta: vector of the ambient temperature. Can be empty [].
%         
%         location: a structure containing the following fields: {'UTC',
%         'latitude','longitude','altitude','n_thetas'}. UTC is a scalar 
%         containing the UTC shift in hours, latitude and longitude are 
%         floats containing GPS coordinates and altitude is float 
%         containing the altitude in meters. n_thetas is a scalar with the
%         desired number of proxies (used only if equallySpaced is false).
%         
%         equallySpaced: if true, the proxies are obtained from a
%         icosahedron mesh, and their orientations are uniformly
%         distributed on a unit sphere. If false, the proxy are computed on
%         a rectangular grid of tilts and azimuts.
%
%         correnction: a string in {'T','TL','TLV'}. 'T' applies a
%         temperature correction, 'TL' applies an additional correction
%         considering th inverter efficiency, 'TLV' consider also the
%         influence of low incidence angles.
% 
% Output: proxy matrix, T*N_p where T is the number of observations in GHI
%         and N_p is the number of proxies 


if nargin<5
    equallySpaced=true;
end

% Define proxies' orientations
if equallySpaced
    [A,T] = generateProxies;
else
    n_thetas_1d = 4; % number of samples over sun azimuth and sun elevation direction. The total number of thetas is (n_thetas_1d)^2
    titl_lims = [5,90]; % from flat to vertical
    az_lims = [90,270]; % from east to west
    [T,A] = meshgrid(linspace(titl_lims(1),titl_lims(2),n_thetas_1d),linspace(az_lims(1),az_lims(2),n_thetas_1d));
end
TA = [T(:),A(:)];
n_proxies=size(TA,1);
[yy,~] = datevec(t);
dy = floor(t-datenum(yy,1,1))+1;
n_obs = size(GHI,1);
proxies = zeros(n_obs,n_proxies);

% Compute the proxies and correct for temperature, inverter efficiency and
% low incidence angles

for i=1:n_proxies
    [proxies(:,i),Ib,Id,Ig,~,~,~,~,~,~,AOI] = irradiance(GHI,t,location.UTC,location.latitude,location.longitude,location.altitude,TA(i,1),TA(i,2),0.1,dy);
    if strcmp(correction,'T')
        proxies(:,i) = correct_proxies(proxies(:,i),Ta);
    elseif strcmp(correction,'TL')
        proxies(:,i) = correct_proxies(proxies(:,i),Ta,1);
    elseif strcmp(correction,'TLV')
        proxies(:,i) = correct_proxies(proxies(:,i),Ta,1,[AOI,Ib,Id,Ig]);
    else
        proxies(:,i) = correct_proxies(proxies(:,i),[]);
    end

end
