function [firings, stim_times, burstiness_trace] = simulate_closed_loop_stim(Ne, Ni, a, b, c, d, S, total_time, stim_neuron_idx, stim_amplitude, stim_duration, burst_thresh)
% Simulates a network of Izhikevich neurons with a closed loop controller
% aiming to turn on stim when a burst is detected
%
% Inputs:
%   ====Network paramters====
%   Ne            - Number of excitatory neurons
%   Ni            - Number of inhibitory neurons
%   a, b, c, d    - Izhikevich parameters (vectors of length Ne+Ni)
%   S             - Synaptic weight matrix (size [Ne+Ni, Ne+Ni])
%   total_time_ms - Duration of simulation in ms
%   ====Stim related parameters=====
%   stim_neuron_idx - affected neuron indices
%   stim_amplitude  - stim amplitude (same for all stim pulses in time)
%   stim_duration   - stim duration in ms (same for all stim pulses in time)
%   ====Control paramters====
%   burst_thresh    - threshold on burstiness signal to send stim 
%
% Output:
%   firings          - [t, neuron_idx] of spike events
%   stim_times       - times of stimulation pulses
%   burstiness_trace - metric used to monitor for bursts and send stim accordingly  

    N = Ne + Ni;
    v = -65 * ones(N,1);
    u = b .* v;
    firings = [];

    % === Parameters ===
    bin_width = 10;  % burstiness window in ms
    step = 1;        % simulation step size in ms
    max_lag = 200;   % for burstiness buffer
    spike_buffer = zeros(max_lag, N);  % rolling window of spikes

    stim_active = false;
    stim_timer = 0;
    stim_times = [];

    % For tracking burstiness value over time
    burstiness_trace = zeros(total_time, 1);

    for t = 1:total_time
        % === Step 1: Compute burstiness ===
        recent_spikes = spike_buffer(max_lag-bin_width+1:end, :);
        pop_spike_count = sum(recent_spikes, 2);
        burstiness = mean(pop_spike_count);  % can change to sum or std

        burstiness_trace(t) = burstiness;

        % === Step 2: Check stim condition ===
        if ~stim_active && burstiness > burst_thresh
            stim_active = true;
            stim_timer = stim_duration;
            stim_times(end+1) = t;
        end

        % === Step 3: Input current ===
        I = [5*randn(Ne,1); 2*randn(Ni,1)];

        if stim_active
            I(stim_neuron_idx) = I(stim_neuron_idx) + stim_amplitude;
            stim_timer = stim_timer - 1;
            if stim_timer <= 0
                stim_active = false;
            end
        end

        % === Step 4: Simulate dynamics ===
        fired = find(v >= 30);
        firings = [firings; t + 0*fired, fired];
        v(fired) = c(fired);
        u(fired) = u(fired) + d(fired);
        I = I + sum(S(:,fired), 2);

        for k = 1:2
            v = v + 0.5 * (0.04*v.^2 + 5*v + 140 - u + I);
        end
        u = u + a .* (b .* v - u);

        % === Step 5: Update spike buffer ===
        spike_row = zeros(1, N);
        spike_row(fired) = 1;
        spike_buffer = [spike_buffer(2:end, :); spike_row];  % shift buffer
    end
end
