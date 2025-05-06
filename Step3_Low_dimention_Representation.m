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
firings = simulate_izhikevich_network(Ne, Ni, a, b, c, d, S, total_time_ms);

%% Calculate PCA and Burstiness
% === PARAMETERS ===
bin_width_ms = 10;
slide_step_ms = 2;
win_ms = 100;               % raster window shown per frame
fps = 20;                   % animation speed
total_time = 2000;          % simulation duration

% === Sliding Window Time Centers ===
t_centers = bin_width_ms/2 : slide_step_ms : total_time - bin_width_ms/2;
n_bins_slide = length(t_centers);

fprintf("Compute Spike Matrix with Sliding Windows\n")
% === Step 1: Compute Spike Matrix with Sliding Windows ===
X_slide = zeros(n_bins_slide, Ne+Ni);  % time bins × neurons

for i = 1:n_bins_slide
    t_start = t_centers(i) - bin_width_ms/2;
    t_end   = t_centers(i) + bin_width_ms/2;
    for j = 1:Ne+Ni
        spikes_j = firings(firings(:,2) == j, 1);
        X_slide(i,j) = sum(spikes_j >= t_start & spikes_j < t_end);
    end
end
fprintf("Calculating PCA\n")
% === Step 2: PCA and Burstiness Score ===
X_slide_z = zscore(X_slide);  % z-score across time
[~, score_slide, ~] = pca(X_slide_z);
X_pca2_slide = score_slide(:,1:2);

pop_count = sum(X_slide, 2);
burstiness_score_slide = zscore(smoothdata(pop_count, 'gaussian', 5));
burstiness_score_slide = burstiness_score_slide(:);  % ensure column

%% Make video 
% Fix colormap range
min_burst = min(burstiness_score_slide);
max_burst = max(burstiness_score_slide);

% === Step 3: Animation ===
figure('Position', [100, 100, 1200, 500]);

% Determine number of frames
step_bins = win_ms / slide_step_ms;  % how many 2ms steps per window
n_frames = n_bins_slide - step_bins;

saveDir = "Step3_Low_dimention_Representation";
if ~exist(saveDir, 'dir')
       mkdir(saveDir)
end

fileName = "LDR";
% Uncomment to save video
% video_writer = VideoWriter(saveDir + filesep + fileName + ".mp4", 'MPEG-4');
% video_writer.FrameRate = fps;
% open(video_writer);

fprintf("starting video \n");
for frame = 1:n_frames
    clf;

    % === Time window for this frame ===
    t_start = t_centers(frame);
    t_end   = t_centers(frame + step_bins);
    
    % === Subplot 1: Spike Raster ===
    subplot(1,2,1);
    spikes_in_win = firings(firings(:,1) >= t_start & firings(:,1) < t_end, :);
    plot(spikes_in_win(:,1), spikes_in_win(:,2), '.k'); hold on;

    % Add shaded green burst intervals
    if exist('burst_intervals', 'var')
        y_limits = ylim;
        for i = 1:size(burst_intervals,1)
            burst_start = max(burst_intervals(i,1), t_start);
            burst_end   = min(burst_intervals(i,2), t_end);
            if burst_end > burst_start
                patch([burst_start burst_end burst_end burst_start], ...
                      [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
                      [0.3 1 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            end
        end
    end
    curr_idx = frame + floor(step_bins/2);
    xline(t_centers(floor(curr_idx)), 'r--')
    xlim([t_start, t_end]);
    ylim([0, Ne+Ni]);
    title(sprintf('Spike Raster (%d–%d ms)', round(t_start), round(t_end)));
    xlabel('Time (ms)'); ylabel('Neuron Index');

    % === Subplot 2: Latent Space ===
    subplot(1,2,2);
    cla;
    scatter(X_pca2_slide(:,1), X_pca2_slide(:,2), 15, [0.8 0.8 0.8], 'filled'); hold on;

    curr_idx = frame + floor(step_bins/2);
    if curr_idx >= 1 && curr_idx <= length(burstiness_score_slide)
        scatter(X_pca2_slide(curr_idx,1), X_pca2_slide(curr_idx,2), ...
                60, burstiness_score_slide(curr_idx), 'filled');
        
    end

    xlabel('PC 1'); ylabel('PC 2');
    title('Latent Space (Colored by Burstiness)');
    colormap(parula);
    colorbar;
    caxis([min_burst, max_burst]);
    axis equal;
    grid on;

    drawnow;

    % Save as video
    % writeVideo(video_writer, getframe(gcf)); % uncomment so save video
end

% close(video_writer); % uncomment so save video

