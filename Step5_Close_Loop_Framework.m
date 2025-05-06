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


total_time = 2000;

stim_neuron_idx = randperm(Ne+Ni, 300); % affected neuron indices
stim_amplitude = -20; % stim pulse amplitude fixed for all pulses
stim_duration = 10; % stim pulse duration fixed for all pulses  
burst_thresh = 5; % threshold for starting the stim

% Uncontrolled model
[firings_unctrl, ~, burstiness_unctrl] = simulate_closed_loop_stim(Ne, Ni, a, b, c, d, S, total_time, [], 0, 0, Inf);

% Controlled model
[firings_ctrl, stim_times_ctrl, burstiness_ctrl] = simulate_closed_loop_stim(Ne, Ni, a, b, c, d, S, total_time, stim_neuron_idx, stim_amplitude, stim_duration, burst_thresh);

%% Making the video
% === Animation parameters ===
win_ms = 500;
step = 10;
t_step = 1;  % ms
n_frames = total_time - win_ms;
saveDir = "Step5_Close_Loop_Framework";
if ~exist(saveDir, 'dir')
       mkdir(saveDir)
end

fileName = "Closeloop system";
% Uncomment to save
% video_writer = VideoWriter(saveDir + filesep + fileName + ".mp4", 'MPEG-4');
% video_writer.FrameRate = fps;
% open(video_writer);

figure('Position', [100, 100, 1200, 600]);

for t = 1:step:n_frames
    clf;

    t_start = t;
    t_end = t + win_ms;

    % --- 1. Uncontrolled Raster with Burst Highlight ---
    subplot(2,2,1);
    spikes = firings_unctrl(firings_unctrl(:,1) >= t_start & firings_unctrl(:,1) < t_end, :);
    plot(spikes(:,1), spikes(:,2), '.k'); hold on;

    % Shade burstiness periods (uncontrolled)
    % Plot burstiness periods in blue in raster (top and bottom)
    highlight_times = t_start:t_end;
    burst_signal_win = burstiness_unctrl(t_start:t_end) > burst_thresh*2;
    burst_time_indices = highlight_times(burst_signal_win);
    
    % Shade each time point where burstiness > threshold
    for b = burst_time_indices
        patch([b b+1 b+1 b], [0 0 Ne+Ni Ne+Ni], ...
              [0.6 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end

    xlim([t_start, t_end]);
    ylim([0, Ne+Ni]);
    title('Uncontrolled Raster with Burst Highlight');

    % --- 2. Uncontrolled Burstiness Trace ---
    subplot(2,2,2);
    plot(t_start:t_end, burstiness_unctrl(t_start:t_end), 'b', 'LineWidth', 2); hold on;
    yline(burst_thresh *2, 'r--', 'Burst Level');
    xlim([t_start, t_end]);
    ylim([min(burstiness_unctrl), max(burstiness_unctrl)]);
    title('Uncontrolled Burstiness');

    % --- 3. Controlled Raster with Burst + Stim Highlight ---
    subplot(2,2,3);
    spikes = firings_ctrl(firings_ctrl(:,1) >= t_start & firings_ctrl(:,1) < t_end, :);
    plot(spikes(:,1), spikes(:,2), '.k'); hold on;

    % Highlight burstiness periods (controlled)
    highlight_times_ctrl = t_start:t_end;
    burst_signal_ctrl = burstiness_ctrl(t_start:t_end) > burst_thresh*2;
    burst_time_indices_ctrl = highlight_times_ctrl(burst_signal_ctrl);
    
    for b = burst_time_indices_ctrl
        patch([b b+1 b+1 b], [0 0 Ne+Ni Ne+Ni], ...
              [0.6 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end

    % % Highlight stim periods
    % for s = 1:length(stim_times_ctrl)
    %     stim_t = stim_times_ctrl(s);
    %     stim_end = stim_t + stim_duration;
    %     if stim_t >= t_start && stim_t <= t_end
    %         patch([stim_t stim_end stim_end stim_t], ...
    %               [0 0 Ne+Ni Ne+Ni], ...
    %               [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    %     end
    % end

    xlim([t_start, t_end]);
    ylim([0, Ne+Ni]);
    title('Controlled Raster with Burst + Burst Highlight');

    % --- 4. Controlled Burstiness Trace ---
    subplot(2,2,4);
    plot(t_start:t_end, burstiness_ctrl(t_start:t_end), 'b', 'LineWidth', 2); hold on;
    yline(burst_thresh * 2, 'r--', 'Burst Level');
    yline(burst_thresh, 'g--', 'Threshold Activating Stim');

    % % Mark stim times
    % stim_in_win = stim_times_ctrl(stim_times_ctrl >= t_start & stim_times_ctrl <= t_end);
    % for s = stim_in_win
    %     xline(s, 'Color', [1 0 0 0.6], 'LineWidth', 1.5);  % transparent red
    % end

    xlim([t_start, t_end]);
    ylim([min(burstiness_ctrl), max(burstiness_ctrl)]);
    title('Controlled Burstiness (Stim marked)');

    drawnow;

    % Save video
    % writeVideo(video_writer, getframe(gcf)); % Uncomment to save
end

% After loop:
% close(video_writer); % Uncomment to save

