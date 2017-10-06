function [DNI,DHI] = engerer2(GHI,Zenit,dayOfYear, t,location,AST)
% This function implements the engerer2 separation model for the DNI DHI
% split. The code is not tested and we are not confident with the results
% this split provides, therefore the disc model is currently used instead.

% ENGERER2 separation model - parameters
c = 4.2336e-2;
b0 = -3.7912;
b1 = 7.5479;
b2 = -1.0036e-2;
b3 = 3.1480e-3;
b4 = -5.3146;
b5 = 1.7073;

HExtra = pvl_extraradiation(dayOfYear);   % extraterrestrial irradiance
% Zenit = deg2rad(90-SunEl);
HExtra_hor = HExtra.*cosd(Zenit);
HExtra_hor_tr = HExtra_hor;
HExtra_hor_tr(Zenit>90)=0;
Kt = GHI./HExtra_hor_tr;

% Compute DKtc
GHI_clear = ghi_clear_sky(t,location,location.UTC); % clear sky global horizontal irradiance
Ktc = GHI_clear./HExtra_hor_tr;

DKtc = Ktc-Kt;
Kde = max(0,1-GHI_clear./GHI);

% Compute Kd and DNI
Kd = c + (1+c)./(1+exp(b0+b1*Kt+b2*AST+b3*Zenit+b4*DKtc)) +b5*Kde;
DHI = Kd.*GHI;
DHI(DHI>GHI)=GHI(DHI>GHI);
DNI = (GHI-DHI)./cosd(Zenit);
DNI(Zenit>90) = 0;
DNI(Zenit>90)=0;

end