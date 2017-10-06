function [SunAz,SunEl] = get_az_el(t,location)
Time = pvl_maketimestruct(t, ones(size(t))*location.UTC); %generate a structured Time
Location = pvl_makelocationstruct(location.latitude,location.longitude,location.altitude); %Generate a structured Location
[SunAz, SunEl, AppSunEl, ~] = pvl_ephemeris(Time,Location); %compute azimuth, elevation (90°-zenit) , apparent elevation and solar time. NB:can take additional arguments
