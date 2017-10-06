function dataset = bin_average(dataset,n_bins)

fields = fieldnames(dataset);
n = length(dataset.(fields{1})); % length of data in dataset

% check length of fields
for i = 1:length(fields)
    if length(dataset.(fields{1})) ~= n
        error('error: All fields must have the same length')
    end
end

n_max = floor(n/n_bins)*n_bins;
for i=1:length(fields)
    var_i = dataset.(fields{i});
    var_i = var_i(1:n_max,:);
    for j=1:size(var_i,2)
        var_i_m(:,j) = mean(reshape(var_i(:,j),n_bins,[]))';
    end
    dataset.(fields{i}) = var_i_m;
    clear var_i_m
end