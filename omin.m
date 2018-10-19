% Minimization function
% Evaluates difference between data and simulation


function oval = omin(p,y01,tspan,tspan2,exp_data)

% Run steady state simulation
[time_ss,species_ss] = ode23s(@oct4f,tspan,y01,[],p);
sz = size(species_ss);
xones = ones(sz);

% Run Diff simulation
y02 = species_ss(241,:);
y02(1)=150;
[time_diff,species_diff] = ode23s(@oct4f,tspan2,y02,[],p);

% Concatenate simulation  results
time_sim = [time_ss;time_diff];
time_sim(241)=[];
species_sim = [species_ss;species_diff];
species_sim(241,:)=[];

% COncatenate data (ones + diff data points)
exp_data=exp_data(:,2:4);
exp_res=[xones(1,1);exp_data'];

time_data = 0:24:72;

% Plot Sim vs Data

plot(time_data,exp_res,'o')
hold
plot(time_sim,species_sim(:,2),'--')
drawnow
hold off 
% 
legend({'Data','Model'})

% Compare steady state and diff simuations against data. Want ONE Value
% represents how good the model simulataed BOTH conditions.

yi = interp1(time_sim,species_sim(:,2),time_data); % get simulated values exactly at data time points

oval = sum((yi- exp_res').^2);

end