function [keep_idx,pv_filtered] = get_clear_day(pv,t,t_filt,rmse_tol,quality_interval)
%Find clear periods using only observations from pv plants. This function uses
%a sliding window approach. Firstly, a low pass filter is applied to the pv
%measurements. Then, the signal are analyzed in a moving window fashion. If
%the mean relative error is less then err_tol in the window, the
%observations are marked as clear observations.
%Inputs:
% pv: a colum vector or matrix of pv signals
% t_filt: filter tau
% sliding_window_size: size of the sliding window. Bigger windows means
% longer selected periods, and fewer periods.
% rel_tol: relative tolerance. The lower, the more selective.

[a,b] = butter(2,1/(t_filt));
pv_filtered = filtfilt(a,b,pv);
pv_filtered(pv_filtered<0)=0;
pv(pv<0)=0;

day_idx = my_days365(t(1),t)+1;
obs_per_day = accumarray(day_idx,ones(size(day_idx,1),1));
obs_per_day = obs_per_day(2);
pv_d = pv(logical(day_idx>2 & day_idx<day_idx(end)));
pv_filtered_d = pv_filtered(logical(day_idx>2 & day_idx<day_idx(end)));

n_max = obs_per_day*floor(size(pv_d,1)/obs_per_day);
pv_d = pv_d(1:n_max);
pv_filtered_d = pv_filtered_d(1:n_max);
pv_d = reshape(pv_d,obs_per_day,[]);
pv_filtered_d = reshape(pv_filtered_d,obs_per_day,[]);

keep_day = mean(((pv_d-pv_filtered_d)./mean(pv_d)).^2).^0.5<rmse_tol & mean(pv_d)>0.5*mean(mean(pv_d(pv_d>0)));
keep_idx = [false(sum(day_idx<=2),1);reshape(repmat(keep_day,obs_per_day,1),[],1);false(sum(day_idx==day_idx(end)),1)];

if ~isempty(quality_interval)
    [yy,mm,dd,hh] = datevec(t(1:length(keep_idx)));
    keep_idx = keep_idx & (quality_interval(1)<hh & hh<quality_interval(2));
end
    
% figure; plot([pv_filtered(keep_idx),pv(keep_idx)])
  