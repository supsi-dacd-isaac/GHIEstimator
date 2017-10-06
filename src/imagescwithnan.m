function [h,hcb] = imagescwithnan(a,cm,nanclr,xdata,ydata,cLims)
% IMAGESC with NaNs assigning a specific color to NaNs

if nargin<6
%# find minimum and maximum
amin=min(a(:));
amax=max(a(:));
else
    amin=cLims(1);
    amax=cLims(2);
end
%# size of colormap
n = size(cm,1);
%# color step
dmap=(amax-amin)/(n-1);

%# standard imagesc
him = imagesc(a,'XData',xdata,'YData',ydata);
%# add nan color to colormap
colormap([nanclr; cm]);
%# changing color limits
caxis([amin-dmap amax]);
%# place a colorbar
hcb = colorbar;
%# change Y limit for colorbar to avoid showing NaN color
ylim(hcb,[amin amax])

if nargout > 0
    h = him;
end