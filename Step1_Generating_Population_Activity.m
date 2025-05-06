clc 
clear all
close all

%% the original parameters from the paper
% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons Inhibitory neurons
Ne=800; Ni=200;
re=rand(Ne,1); ri=rand(Ni,1);
a=[0.02*ones(Ne,1); 0.02+0.08*ri];
b=[0.2*ones(Ne,1); 0.25-0.05*ri];
c=[-65+15*re.^2; -65*ones(Ni,1)];
d=[8-6*re.^2; 2*ones(Ni,1)];
S=[0.5*rand(Ne+Ni,Ne), -rand(Ne+Ni,Ni)];
v=-65*ones(Ne+Ni,1); % Initial values of v
u=b.*v; % Initial values of u

total_time_ms = 2000;
firing = simulate_izhikevich_network(Ne, Ni, a, b, c, d, S, total_time_ms);
experiment = "Original Parameters";
title(experiment)

saveDir = "Step1_Generating_Population_Activity";
if ~exist(saveDir, 'dir')
       mkdir(saveDir)
end

fileName = experiment;
% savefig(saveDir + filesep + fileName)               % uncomment to save
% saveas(gcf, saveDir + filesep + fileName+".png")    % uncomment to save


%% Increased excitatory neurons 
Ne=1000; Ni=200;
re=rand(Ne,1); ri=rand(Ni,1);
a=[0.02*ones(Ne,1); 0.02+0.08*ri];
b=[0.2*ones(Ne,1); 0.25-0.05*ri];
c=[-65+15*re.^2; -65*ones(Ni,1)];
d=[8-6*re.^2; 2*ones(Ni,1)];
S=[0.5*rand(Ne+Ni,Ne), -rand(Ne+Ni,Ni)];
v=-65*ones(Ne+Ni,1); % Initial values of v
u=b.*v; % Initial values of u

total_time_ms = 2000;
firing = simulate_izhikevich_network(Ne, Ni, a, b, c, d, S, total_time_ms);
title("Increased excitatory neurons")
experiment = "Increased excitatory neurons e=" + string(Ne)...
          + " i=" + string(Ni);
title(experiment)

saveDir = "Step1_Generating_Population_Activity";
if ~exist(saveDir, 'dir')
       mkdir(saveDir)
end

fileName = experiment;
% savefig(saveDir + filesep + fileName)               % uncomment to save
% saveas(gcf, saveDir + filesep + fileName+".png")    % uncomment to save

%% Increased inhibitory neurons 
Ne=800; Ni=400;
re=rand(Ne,1); ri=rand(Ni,1);
a=[0.02*ones(Ne,1); 0.02+0.08*ri];
b=[0.2*ones(Ne,1); 0.25-0.05*ri];
c=[-65+15*re.^2; -65*ones(Ni,1)];
d=[8-6*re.^2; 2*ones(Ni,1)];
S=[0.5*rand(Ne+Ni,Ne), -rand(Ne+Ni,Ni)];
v=-65*ones(Ne+Ni,1); % Initial values of v
u=b.*v; % Initial values of u

total_time_ms = 2000;
firing = simulate_izhikevich_network(Ne, Ni, a, b, c, d, S, total_time_ms);
title("Increased inhibitory neurons")
experiment = "Increased inhibitory neurons e=" + string(Ne)...
          + " i=" + string(Ni);
title(experiment)

saveDir = "Step1_Generating_Population_Activity";
if ~exist(saveDir, 'dir')
       mkdir(saveDir)
end

fileName = experiment;
% savefig(saveDir + filesep + fileName)               % uncomment to save
% saveas(gcf, saveDir + filesep + fileName+".png")    % uncomment to save

%% Add certaion types of neurons such as chattering neurons by modulating 
% model parameters 
Ne=850;                 Ni=350;
re=rand(Ne,1);          ri=rand(Ni,1);
a=[0.02*ones(Ne,1);     0.02+0.08*ri];
b=[0.2*ones(Ne,1);      0.25-0.05*ri];
c=[-55+20*re.^2;        -55*ones(Ni,1)];
d=[8-6*re.^2;           2*ones(Ni,1)];
S=[0.45*rand(Ne+Ni,Ne),  -1.1*rand(Ne+Ni,Ni)]; 
v=-65*ones(Ne+Ni,1);    % Initial values of v
u=b.*v;                 % Initial values of u

total_time_ms = 2000;
firing = simulate_izhikevich_network(Ne, Ni, a, b, c, d, S, total_time_ms);
title("BG like population (chattering)")
experiment = "BG like population";
title(experiment)

saveDir = "Step1_Generating_Population_Activity";
if ~exist(saveDir, 'dir')
       mkdir(saveDir)
end

fileName = experiment;
% savefig(saveDir + filesep + fileName)               % uncomment to save
% saveas(gcf, saveDir + filesep + fileName+".png")    % uncomment to save

