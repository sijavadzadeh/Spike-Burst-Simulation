function [firings, stim_log] = simulate_izhikevich_with_stimulation(Ne, Ni, a, b, c, d, S, total_time_ms, stim)
% Simulates Izhikevich network with multiple stim events passed as arrays
% Inputs:
%   ====Network paramters====
%   Ne            - Number of excitatory neurons
%   Ni            - Number of inhibitory neurons
%   a, b, c, d    - Izhikevich parameters (vectors of length Ne+Ni)
%   S             - Synaptic weight matrix (size [Ne+Ni, Ne+Ni])
%   total_time_ms - Duration of simulation in ms
%   ====Stim related parameters=====
%   stim must contain: 
%   .stim_start       - [1 x n] vector of stim start times
%   .stim_duration    - [1 x n] vector of stim durations
%   .stim_amplitude   - [1 x n] vector of stim amplitudes
%   .stim_neuron_idx  - [1 x m] vector of neuron indices affected by stim
%
% Output:
%   firings          - [t, neuron_idx] of spike events
%   stim_log         - logs stimulation periods


    N = Ne + Ni;
    v = -65 * ones(N,1);
    u = b .* v;
    firings = [];

    stim_log = zeros(total_time_ms, N);

    n_stim_events = length(stim.stim_start);

    for t = 1:total_time_ms
        I = [5*randn(Ne,1); 2*randn(Ni,1)];
        
        % Modulate I based on stim parameters 
        for s = 1:n_stim_events
            stim_start = stim.stim_start(s);
            stim_end   = stim_start + stim.stim_duration(s) - 1;

            if t >= stim_start && t <= stim_end
                idx = stim.stim_neuron_idx;
                I(idx) = I(idx) + stim.stim_amplitude(s);
                stim_log(t, idx) = stim_log(t, idx) + stim.stim_amplitude(s);
            end
        end

        fired = find(v >= 30);
        firings = [firings; t + 0*fired, fired];

        v(fired) = c(fired);
        u(fired) = u(fired) + d(fired);
        I = I + sum(S(:,fired), 2);

        for step = 1:2
            v = v + 0.5 * (0.04*v.^2 + 5*v + 140 - u + I);
        end
        u = u + a .* (b .* v - u);
    end
end
