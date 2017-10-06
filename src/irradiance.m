function [I,Ib,Id,Ig,DNI,DHI,AM,HExtra,SunZen,SunAz,AOI]=irradiance(GHI,t,UTC,L_lat,L_long,L_alt,L_tilt,L_az,L_alb,dayOfYear)
%% this function computes the solar irradiance over an oriented surface, given the following inputs:
%
% -GHI global horizontal irradiance
% -t vector of time in datenum format, realtive to GHI measurements
% -UTC offset with respect to Greenwich time; between 0-23
% -L_lat location latitude
% -L_long location longitude
% -L_alt location altitude
% -Ltilt location tilt for the oriented surface
% -L_az location azimuth for the oriented surface
% -L_alb location albedo
%
% The function computes the total irradiance collected by an oriented surface
% Total irradiance I is calculated as:
% I=Ib+Id+Ig
% where Ib is the direct beam irradiance, which is Eb=DNI*cos(AOI)
% Id is the diffuse irradiance on the oriented surface, computed through perez model
% Ig is the ground contribution, computed as Ig=GHI*Albedo*(1-cos(SurfTilt))/2
%
%  Outputs:
% -I total irradiance
% -Ib direct beam irradiance
% -Id diffuse irradiance
% -Ig irradiance from the ground
% -DNI direct normal irradiance
% -DHI diffuse horizontal irradiance
% -AM air mass
% -HExtra extra terrestrial irradiance
% -SunZen zenit angle of the sun
% -SunAz  azimuth angle of the sun


%% Split GHI in DNI and DHI with disc model

Time = pvl_maketimestruct(t, ones(size(t))*UTC);                % generate a structured Time
Location = pvl_makelocationstruct(L_lat,L_long,L_alt);          % generate a structured Location
[SunAz, SunEl, AppSunEl, AST] = pvl_ephemeris(Time,Location);   % compute azimuth, elevation (90-zenit) , apparent elevation and solar time. NB:can take additional arguments
Zenit=90-SunEl;

% [DNI,DHI] = engerer2(GHI,Zenit,dayOfYear,t,location,AST);

DNI = pvl_disc(GHI,Zenit, dayOfYear);                           % DNI from disc model NB:can take additional arguments
DHI = GHI - cosd(Zenit).*DNI;                                   % DHI form disc model

%% Compute projection with Perez model

surfAz=ones(length(Time.second),1).*L_az;                       % generate vectors for parameters
surfTilt=ones(length(Time.second),1).*L_tilt;                   % generate vectors for parameters
Albedo=ones(length(Time.second),1).*L_alb;                      % generate vectors for parameters
HExtra = pvl_extraradiation(dayOfYear);                         % extraterrestrial irradiance
SunZen = 90 - AppSunEl;                                         % real zenith
AM = pvl_relativeairmass(SunZen);                               % air mass
AM(isnan(AM)) = 20;
Id = pvl_haydavies1980(surfTilt, surfAz, DHI, DNI, HExtra, SunZen, SunAz);

%% Computation of Ig and Ib

Ig= pvl_grounddiffuse(surfTilt, GHI, Albedo);                   % irradiance from ground
AOI = pvl_getaoi(surfTilt, surfAz, SunZen, SunAz);              % get AOI
Ib = 0*AOI;                                                     % initiallize variable
Ib(AOI<90) = DNI(AOI<90).*cosd(AOI(AOI<90));                    % calculate only when sun is in view of the plane of array

%% Compute total Irradiation
I=Id+Ib+Ig;


