function out = plot_cols(i)
    col = [0    0.4470    0.7410;  % Blue
        0.8500    0.3250    0.0980; % Red
        0.9290    0.6940    0.1250; % Yellow
        0.4940    0.1840    0.5560; % Purple
        0.4660    0.6740    0.1880; % Green
        0.3010    0.7450    0.9330; % Light Blue
        0.6350    0.0780    0.1840];% Dark Red
    
    i = mod(i-1,length(col(:,1)))+1;
    out = col(i,:);
end