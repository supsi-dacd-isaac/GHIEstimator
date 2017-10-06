function [az,tilt,el]=generateProxies
% Generate a base icosahedron mesh
TR=IcosahedronMesh;
 
% Subvivide the base mesh
for i=2
    TR=SubdivideSphericalMesh(TR,1);
    [azimuth,elevation]=cart2sph(TR.X(:,1),TR.X(:,2),TR.X(:,3));
    azdeg=azimuth/pi*180+180;
    eldeg=elevation/pi*180;
    %exclude north (this is so far hardcoded)
    sel=(azdeg>=50&azdeg<=310&eldeg>=0)| eldeg>50;
    sum(sel);
    az=azdeg(sel);
    el=eldeg(sel);
    tilt=90-el;
end