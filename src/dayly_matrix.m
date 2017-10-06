function [d_mat,t,obs_per_day] = dayly_matrix(x,t_datenum,quality_interval)
%take the x signal associated with the t_datenum vect
%and return d_mat of size n_days*n_obs_per_day. Note that n_days*n_obs_per_day
% is likely to be different from length(x).

if size(x,2)>size(x,1)
    fprintf('Matrix must be slim')
else
    day_idx = my_days365(t_datenum(1),t_datenum)+1;
    [~,~,~,hh] = datevec(t_datenum);
    keep_these = logical(day_idx>2 & day_idx<day_idx(end));
    if  ~isempty(quality_interval)
        keep_these = keep_these & (quality_interval(1)<hh & hh<quality_interval(2));
    end
    
    obs_per_day = unique(accumarray(day_idx,keep_these));
    
    if  length(obs_per_day)==1
        obs_per_day = obs_per_day(1); 
    else 
        obs_per_day = obs_per_day(2);
    end
    
    
    if size(x,2)==1
        x = x(keep_these);
        d_mat = reshape(x,obs_per_day,[]);
    else
        for i=1:size(x,2)
            x_i = x(keep_these,i);
            d_mat(:,:,i) = reshape(x_i,[],obs_per_day);
        end
    end
    t = t_datenum(keep_these);
end