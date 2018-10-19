% get OCT4 time course data
time_data = W; % time points
species_data = Z; % species

% get qtpcr data
file_name = 'oct4_cdx2_timecourse';
[num,txt,raw] = xlsread(file_name,'OCT4 - ACTB');

oct4_data = num(:,12);
octavg = mean(oct4_data(1:3));
norm_oct4=oct4_data/octavg;
time_zero= mean(norm_oct4(1:3));
time_one = mean(norm_oct4(4:6));
time_two = mean(norm_oct4(7:9));
time_three = mean(norm_oct4(10:12));
exp_data = [time_zero, time_one, time_two, time_three];

% Plot data
figure
h = gca;
h.XLabel.String = 'Time (min)';
h.YLabel.String = 'Fold Change';

%initial conditions. 
diff_signal = 1;
x1 = 1; %oct4
x2 = 1;
y1 = 1; %cdx2
y2 = 1;
% smax
% s1 = 0.6914;
% s2 = 2.0906;
s1 = .9;
s2 = .9;
%concentration at which the response is at half max
K1 = 1.1;
K2 = 1.1;

%Hill coefficient

n1 = 2.2;
n2 = 2.2;

% protein synthesis constant
beta1 = 1.2206;
beta2 = 1.2273;

%degradation constant.
a1 = 1.27;
a3 = 1.3;
a4 = 1.03;

tspan1 = 0:0.1:24;
tspan2 = 24:0.1:72;

A = [];
b = [];
Aeq = [];
beq = [];
% lb = [.001,.001,.001,.001,.001,.001,.001,.001,log(2)/20,log(2)/20,log(2)/20];
lb = log10([.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001]);
ub = log10([100,100,100,100,100,100,100,100,100,100,100]);
% ub = [100,100,100,100,100,100,100,100,log(2)/5,log(2)/5,log(2)/5];

y01 = [diff_signal,x1,x2,y1,y2]; %initial conditions
y02 = [150,x1,x2,y1,y2]; %initial conditions

p0 = log10([s1,s2,K1,K2,n1,n2,beta1,beta2,a1,a3,a4]); % initial guesses for parameters
options = optimset('Display','iter-detailed');

p = fminsearch(@omin,p0,options,y01,tspan1,...
                                 tspan2,exp_data);
warning off 


