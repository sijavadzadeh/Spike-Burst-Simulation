function firings = simulate_izhikevich_network(Ne, Ni, a, b, c, d, S, total_time_ms)
% Simulates a network of Izhikevich neurons
%
% Inputs:
%   Ne            - Number of excitatory neurons
%   Ni            - Number of inhibitory neurons
%   a, b, c, d    - Izhikevich parameters (vectors of length Ne+Ni)
%   S             - Synaptic weight matrix (size [Ne+Ni, Ne+Ni])
%   total_time_ms - Duration of simulation in ms
%
% Output:
%   firings       - [t, neuron_idx] of spike events

    if nargin < 7
        total_time_ms = 5000;  % default duration
    end

    N = Ne + Ni;
    v = -65 * ones(N,1);   % Membrane potential
    u = b .* v;            % Recovery variable
    firings = [];          % To store spike events

    for t = 1:total_time_ms
        I=[5*randn(Ne,1);2*randn(Ni,1)]; % thalamic input
        fired=find(v>=30); % indices of spikes
        firings=[firings; t+0*fired,fired];
        v(fired)=c(fired);
        u(fired)=u(fired)+d(fired);
        I=I+sum(S(:,fired),2);
        v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
        v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
        u=u+a.*(b.*v-u); % stability
    end
    
    % Plot spiking activity 
    figure('Position', [100, 100, 1200, 500])
    plot(firings(:,1),firings(:,2),'.')
    ylabel("Neuron #")
    xlabel("Time(ms)")
    yline(Ne, 'r-')
    yticks([1, Ne + Ni])
    xl = xlim;
    box off
     
    % Inhibitory neurons at top 
    x_pos = xl(1) - (xl(2)-xl(1))*0.04; 
    text(x_pos, Ne, 'Inhibitory', 'Color', [0.6 0.1 0.1], ...
         'FontSize', 9, 'HorizontalAlignment', 'left','Rotation', 90);
    % Excitatory neurons at bottom
    text(x_pos, Ne/2, 'Excitatory', 'Color', [0.1 0.1 0.6], ...
         'FontSize', 9, 'HorizontalAlignment', 'right','Rotation', 90);
end