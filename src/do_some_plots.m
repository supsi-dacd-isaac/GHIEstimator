function do_some_plots(GHI,GHI_est,t_filt,GHI_sat,GHI_sat2,saveplot,title_str)

if nargin<6
    title_str = '';
end

% retrieve errors
est_err_ghi = GHI-GHI_est;

sampling_time = mean(diff(t_filt(1:10)))*3600*24;

% Time series of GHI for the first d_select days
d_select = 3;
sampling_per_day = round(3600*24/sampling_time);
select = [1:1:min(d_select*sampling_per_day,length(GHI))];

figure('position',[0,0,1000,800]);
plot(t_filt(select),[GHI(select),GHI_est(select)]);hold on;
legend_labs = {'GHI_{obs}','GHI_{est}'};

% Plot satellite based GHI signals if present in the dataset
if ~isempty(GHI_sat)
    plot(t_filt(select),GHI_sat(select),'linewidth',2);
    legend_labs{end + 1} = 'GHI_{sat}';
end
if ~isempty(GHI_sat2)
    plot(t_filt(select),GHI_sat2(select),'linewidth',2);
    legend_labs{end + 1} = 'GHI_{sat 2}';
end
    
legend(legend_labs)
xlabel('GHI [W/m2]')
datetick
title(strcat(title_str,' Time series'))
if saveplot
    saveas(gcf,'fig/TimeSeries.png')
end


% error histograms
n_bins = 1000;
figure('position',[0,0,1000,800]);
histogram(est_err_ghi,n_bins,'normalization','pdf');hold on;

legend_labs = {'Estimated'};

if ~isempty(GHI_sat)
sat_err_ghi = GHI-GHI_sat;
histogram(sat_err_ghi,n_bins,'normalization','pdf');
legend_labs{end+1} = 'Satellite';
end

if ~isempty(GHI_sat2)
sat_err_ghi = GHI-GHI_sat2;
histogram(sat_err_ghi,n_bins,'normalization','pdf');
legend_labs{end+1} = 'Satellite 2';
end

legend(legend_labs)
title(strcat(title_str,' Err dist'))
xlim([-200,200])
ylabel('pdf')
if saveplot
    saveas(gcf,'fig/ERR.png')
end

