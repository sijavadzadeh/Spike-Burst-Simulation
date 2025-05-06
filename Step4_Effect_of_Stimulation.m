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

total_time_ms = 2000;

% Define stim struct with arrays
stim.stim_start = randperm(total_time_ms, 100); % stim pulse start times
stim.stim_duration = ones(1, length(stim.stim_start)) + 20; % stim pulse duration for each pulse
stim.stim_amplitude = ones(1, length(stim.stim_start)) + 20; % stim pulse amplitude for each pulse
stim.stim_neuron_idx = randperm(Ne+Ni, 3);  % affected neurons 

% Run simulation
[firings, stim_log] = simulate_izhikevich_with_stimulation(Ne, Ni, a, b, c, d, S, total_time_ms, stim);

%% === Basic raster plot ===
figure;
plot(firings(:,1), firings(:,2), '.k'); hold on;

xlabel('Time (ms)');
ylabel('Neuron Index');
title('Raster Plot with Stimulation 50 neurons');

% === Add shaded rectangles for each stimulation event ===
y_limits = ylim;  % use full vertical range of neuron indices

for i = 1:length(stim.stim_start)
    stim_start = stim.stim_start(i);
    stim_end = stim_start + stim.stim_duration(i);

    % Draw transparent rectangle for this stim
    patch([stim_start stim_end stim_end stim_start], ...
          [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
          [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);  % yellow-orange
end


%% PSTH for affected neurons
% === PARAMETERS ===
pre_time = 200;   % ms before stim
post_time = 400;  % ms after stim
t_window = -pre_time : post_time;  % relative time axis

% For each neuron that was stimulated in at least one stim event
all_stimulated_neurons = stim.stim_neuron_idx;
n_neurons_to_plot = length(all_stimulated_neurons);

for n = 1:n_neurons_to_plot
    neuron_id = all_stimulated_neurons(n);
    stim_times_for_this_neuron = [];

    % Collect all stim onsets that targeted this neuron
    for s = 1:length(stim.stim_start)
        if ismember(neuron_id, stim.stim_neuron_idx)
            stim_times_for_this_neuron(end+1) = stim.stim_start(s);
        end
    end

    % Create raster data aligned to those stim times
    figure('units', 'normalized','outerposition',[0.4 0.4 0.55 0.45])
    hold on;
    title(sprintf('Neuron %d - Stim-triggered Raster', neuron_id));
    xlabel('Time from stim onset (ms)');
    ylabel('Stimulation repetition');
    
    for trial = 1:length(stim_times_for_this_neuron)
        stim_time = stim_times_for_this_neuron(trial);

        % Extract spikes for this neuron around stim time
        t_min = stim_time - pre_time;
        t_max = stim_time + post_time;

        spikes = firings(firings(:,2) == neuron_id, 1);
        aligned_spikes = spikes(spikes >= t_min & spikes <= t_max) - stim_time;

        % Plot row for this trial
        y_val = trial;
        for spike_time = aligned_spikes'
            plot(spike_time, y_val, '.k');
        end
    end

    xlim([-pre_time, post_time]);
    ylim([0, length(stim_times_for_this_neuron)+1]);
    hold on
    plot([0,stim.stim_duration(1)], [0,0], "Color", "red", LineWidth=5)

    saveDir = "Step4_Effect_of_Stimulation";
    if ~exist(saveDir, 'dir')
           mkdir(saveDir)
    end
    
    fileName = "Pos Amp PSTH neuron " + string(neuron_id);
    % savefig(saveDir + filesep + fileName)             % uncomment to save
    % saveas(gcf, saveDir + filesep + fileName+".png")  % uncomment to save
end


%% Calculate low dimensional representation based on all firing rates
% === Sliding Window Parameters ===
bin_width_ms = 10;
step_ms = 2;
total_time = max(firings(:,1));  % or use simulation duration

% Time centers for bins
t_centers = bin_width_ms/2 : step_ms : total_time - bin_width_ms/2;
n_bins_slide = length(t_centers);
n_neurons = Ne + Ni;

% === Step 1: Sliding spike count matrix ===
X_slide = zeros(n_bins_slide, n_neurons);  % rows = time bins, cols = neurons

for i = 1:n_bins_slide
    t_start = t_centers(i) - bin_width_ms/2;
    t_end   = t_centers(i) + bin_width_ms/2;

    for j = 1:n_neurons
        spikes_j = firings(firings(:,2) == j, 1);
        X_slide(i,j) = sum(spikes_j >= t_start & spikes_j < t_end);
    end
end

% === Step 2: Normalize ===
X_slide_z = zscore(X_slide);  % z-score across time (row-wise)

% === Step 3: PCA ===
[~, score, ~] = pca(X_slide_z);  % score: rows = time, cols = PC
X_pca2_slide = score(:,1:2);     % 2D latent trajectory

%% create a flow field that shows effect of stim in LDR
% Parameters
arrow_step = 5;  % skip every N points to avoid clutter

% Compute flow vectors (differences)
X = X_pca2_slide(:,1);
Y = X_pca2_slide(:,2);
dX = [diff(X); 0];
dY = [diff(Y); 0];

% Plot latent space with flow arrows
figure;
hold on;

% Plot background trajectory
plot(X, Y, '-', 'Color', [0.8 0.8 0.8]);

% Color-coded arrows based on stim timing
for i = 1:arrow_step:length(X)-1
    t_now = t_centers(i);

    % Determine state: pre, stim, or post
    is_stim = false;
    for s = 1:length(stim.stim_start)
        stim_start = stim.stim_start(s);
        stim_end   = stim_start + stim.stim_duration(s);
        if t_now >= stim_start && t_now <= stim_end
            is_stim = true;
            break;
        end
    end

    if is_stim
        arrow_color = [1 0.2 0.2];  % red for stim
        quiver(X(i), Y(i), dX(i), dY(i), 0, 'Color', arrow_color, 'LineWidth', 1);
    end

    
end

xlabel('PC 1'); ylabel('PC 2');
title('Latent State Flow with Stimulation');
axis equal;
grid on;
legend({'Trajectory', 'Flow arrows (stim effect)'});



