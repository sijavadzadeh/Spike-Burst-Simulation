clc
clear all

%% Simulate spiking activity (BG-like set of paramters from step 1)
Ne=800;                 Ni=350;
re=rand(Ne,1);          ri=rand(Ni,1);
a=[0.02*ones(Ne,1);     0.02+0.08*ri];
b=[0.2*ones(Ne,1);      0.25-0.05*ri];
c=[-55+20*re.^2;        -55*ones(Ni,1)];
d=[8-6*re.^2;           2*ones(Ni,1)];
S=[0.45*rand(Ne+Ni,Ne),  -1.1*rand(Ne+Ni,Ni)];
v=-65*ones(Ne+Ni,1);    % Initial values of v
u=b.*v;                 % Initial values of u

total_time_ms = 5000;
firings = simulate_izhikevich_network(Ne, Ni, a, b, c, d, S, total_time_ms);


%% Burst detection using spike counts in the population  
% ------------ PARAMETERS ------------
bin_size = 5;           % ms
total_time = 5000;       % total simulation time in ms
std_thresh = 1;          % number of standard deviations above mean
min_consec_bins = 2;     % minimum consecutive high bins to define a burst
edges = 0:bin_size:total_time;

% ------------ BIN SPIKES ------------
bin_counts = histcounts(firings(:,1), edges);
n_bins = length(bin_counts);

% ------------ COMPUTE BURST THRESHOLD ------------
% mu = mean(bin_counts);
% sigma = std(bin_counts);
mu = median(bin_counts); % More robust estimates based on median 
sigma = median(abs(bin_counts-mu)) * 1.4826; % More robust estimates based on median 
threshold = mu + std_thresh * sigma;

% ------------ INITIALIZE OUTPUTS ------------
burst_binary = false(n_bins, 1);  % 1 if bin is part of a burst
burst_intervals = [];            % stores [start_time, end_time] per burst


% ------------ BURST DETECTION (Robust) ------------
high_bins = find(bin_counts > threshold);
burst_binary = false(n_bins, 1);
burst_intervals = [];

if ~isempty(high_bins)
    runs = {};   % cell array to collect runs
    current_run = high_bins(1);

    for i = 2:length(high_bins)
        if high_bins(i) == high_bins(i-1) + 1
            current_run = [current_run, high_bins(i)];
        else
            runs{end+1} = current_run;
            current_run = high_bins(i);
        end
    end
    runs{end+1} = current_run;  % add final run
    
    % Filter runs and extract binary + intervals
    for i = 1:length(runs)
        if length(runs{i}) >= min_consec_bins
            idx_range = runs{i};
            burst_binary(idx_range) = true;
            burst_start_time = edges(idx_range(1));
            burst_end_time = edges(idx_range(end) + 1);  % end of last bin
            burst_intervals = [burst_intervals; burst_start_time, burst_end_time];
        end
    end
end


% ------------ VISUALIZATION ------------
figure('units', 'normalized','outerposition',[0 0 1 1])
axx(1) = subplot(2,1,1);
plot(firings(:,1), firings(:,2), '.k');
hold on;
for i = 1:size(burst_intervals,1)
    xline(burst_intervals(i,1), 'r--', 'LineWidth', 1.2);
end
title('Raster Plot with Burst Onsets');
xlabel('Time (ms)');
ylabel('Neuron Index');
yline(Ne, 'r-')
yticks([1, Ne + Ni])
xl = xlim;
box off
 
% Inhibitory neurons at top 
x_pos = xl(1) - (xl(2)-xl(1))*0.01; 
text(x_pos, Ne, 'Inhibitory', 'Color', [0.6 0.1 0.1], ...
     'FontSize', 12, 'HorizontalAlignment', 'left','Rotation', 90);
% Excitatory neurons at bottom
text(x_pos, Ne/2, 'Excitatory', 'Color', [0.1 0.1 0.6], ...
     'FontSize', 12, 'HorizontalAlignment', 'right','Rotation', 90);

y_limits = ylim;  % Get current y-axis range

for i = 1:size(burst_intervals,1)
    x_start = burst_intervals(i,1);
    x_end   = burst_intervals(i,2);
    width   = x_end - x_start;

    % Draw rectangle behind the spikes
    patch('XData', [x_start x_end x_end x_start], ...
          'YData', [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
          'FaceColor', [0.3 1 0.3], ...      % light green
          'EdgeColor', 'none', ...
          'FaceAlpha', 0.2);                 % transparency
end


axx(2) = subplot(2,1,2);
bar(edges(1:end-1), bin_counts, 'k');
hold on;
yline(threshold, 'r--', 'LineWidth', 1.5);
for i = 1:size(burst_intervals,1)
    xline(burst_intervals(i,1), 'r--');
end
xlabel('Time (ms)');
ylabel('Spike Count / Bin');
title(sprintf('Population Firing Rate (Threshold = Mean + %d*STD)', std_thresh));
linkaxes(axx,'x')


saveDir = "Step2_Detecting_Bursts";
if ~exist(saveDir, 'dir')
       mkdir(saveDir)
end

fileName = "Detecting_Bursts";
% savefig(saveDir + filesep + fileName)             % uncomment to save
% saveas(gcf, saveDir + filesep + fileName+".png")  % uncomment to save

