function d = my_days365(d1,d2) 
%DAYS365 Days between dates based on 365 day year. 
%   D = DAYS365(D1,D2) returns the number of days between D1 and D2 based 
%   on a 365 day year. .  Enter dates as serial date numbers date strings
%   or datetime arrays.  If D2 is less than D1, the value returned will be 
%   negative.  
% 
%   For example, d = days365('1-Dec-1993', '12-Jan-1994') 
%   or d = days365(728264,728306) returns d = 42. 
% 
%   See also DAYS360, DAYSACT, DAYSDIF, DATENUM.

%       Copyright 1995-2015 The MathWorks, Inc.
 
if nargin < 2 
  error(message('finance:days365:missingInputs')) 
end
if ischar(d1) || ischar(d2) || isdatetime(d1) || isdatetime(d2)
  d1 = datenum(d1); 
  d2 = datenum(d2); 
end 
sz = [size(d1);size(d2)]; 
if length(d1) == 1 
  d1 = d1*ones(max(sz(:,1)),max(sz(:,2))); 
end 
if length(d2) == 1 
  d2 = d2*ones(max(sz(:,1)),max(sz(:,2))); 
end 
sizes = [size(d1);size(d2)];
if any(sizes(:,1)~=sizes(1,1)) || any(sizes(:,2)~=sizes(2,2))
  error(message('finance:days365:invalidInputDims'))
end
 
% Need cumulative sum of days at beginning of each month 
% to convert months into days. 
daytotal = [0;31;59;90;120;151;181;212;243;273;304;334]; 
 
c1 = datevec(d1(:)); 
c2 = datevec(d2(:)); 
 
tempd = 365 * (c2(:, 1) - c1(:, 1)) + daytotal(c2(:, 2)) -... 
            daytotal(c1(:, 2)) + c2(:, 3) - c1(:, 3); 
d = reshape(tempd,size(d1));
