function long_sequences_idx = find_long_sequences(x,sliding_window_size)
% Find sequences of sequences longher than sliding_window_size and remove them
% x is a matrix of temporal signals, size t x n_var

ix=[];
long_sequences_idx=false(size(x));
for p=1:length(x)
    if x(p)
        ix=[ix p];
    else
        if length(ix)>=sliding_window_size
            long_sequences_idx(ix)=true;
        end
        ix=[];
    end
end