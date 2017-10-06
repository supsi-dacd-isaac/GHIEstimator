function GHI_clearsky = ghi_clear_sky(t,location,UTC)

Time = pvl_maketimestruct(t, ones(size(t))*UTC); %generate a structured Time
Location = pvl_makelocationstruct(location.latitude,location.longitude,location.altitude); %Generate a structured Location
[SunAz, SunEl, AppSunEl, ~] = pvl_ephemeris(Time,Location); %compute azimuth, elevation (90�-zenit) , apparent elevation and solar time. NB:can take additional arguments
ApparentZenith=90-AppSunEl;
GHI_clearsky = pvl_clearsky_haurwitz(ApparentZenith);