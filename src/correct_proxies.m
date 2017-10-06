function proxy = correct_proxies(proxy,T,lowG,CorrV)
% Correct the proxy signals considering temperature, inverter efficiency
% and low incidence angles. 
% Inputs: proxy: a T*N_p matrix, where T is the total number of
%                observations and N_p is the number of proxies. 
%         T: vector of ambient temperature. Can be emptsy [].
%         lowG: boolean, if true a correction for the inverter efficiency
%         is applied.
%         corrV: a T*4  matrix, which columns contains in this order: angle
%         of incidence, beam component of irradiance, diffuse component of
%         irradiance and ground reflected component of irradiance.
% 
% Output: proxy: matrix of corrected proxies.

if nargin<3
    lowG=0;
end

if nargin<4
    CorrV=[];
end

if ~isempty(CorrV)
    b0=0.05;
    AOI = CorrV(:,1);
    Ib = CorrV(:,2);
    Id = CorrV(:,3);
    Ig = CorrV(:,4);
    IAM = max(1  - b0*(1./cos(min(AOI(:,1)*pi/180,pi/2)) - 1),0);
    proxy = IAM.*Ib +0.95*(Id+Ig); 
end

% Correct the Ppv proxy considering the effect of T
if ~isempty(T)
    gamma = -0.34/100;
    T_nom = 25;
    T_bom = T+proxy*(0.0284+3e-3);
    
    % correct the proxy
    proxy = proxy.*(1+gamma*(T_bom-T_nom));
end

if lowG
    k2 = 0.942; 
    k3 = -5.02e-2; 
    k4 = -3.77e-2;
    Istc = 1000;
    
    % fitted from typical inverter and polycrystalline module data
    Eff = max(0.1, k2 + k3*log(proxy/Istc) + k4*log(proxy/Istc).^2);
    
    proxy = proxy.*Eff;
end

