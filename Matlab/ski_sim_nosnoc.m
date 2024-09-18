% Ski jump system simulation

clear all;
clc;
close all

% REQUIRES CASADI AND NOSNOC IN MATLAB PATH
% Publication results obtained with: casadi-v3.5.5, nosnoc_v050
addpath('YOURPATH\casadi-windows-matlabR2016a-v3.5.5')
addpath('YOURPATH\nosnoc-main_v050\src')

import casadi.*

%% Flags
b_ext = 0;

%%
N_FE = 3;
N_stages = 30;
T_sim = 5;
%% Initialize nosnoc
problem_options = NosnocProblemOptions();
solver_options = NosnocSolverOptions();

%% Settings
problem_options.time_freezing = 1;
problem_options.time_freezing_inelastic = 1;
problem_options.stagewise_clock_constraint = 0;
problem_options.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
solver_options.print_level = 3;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_stages;
solver_options.opts_casadi_nlp.ipopt.max_iter = 3e3;
% settings.homotopy_update_rule = 'superlinear';

%% Initialize model
[model, plot_fc1, plot_fc2] = ski_model_sim_nosnoc();

%% Solve
integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();

%% Unfold results
unfold_struct(results,'base');
if ~b_ext
    q1 = results.x(1,:);
    q2 = results.x(2,:);
    dq1 = results.x(3,:);
    dq2 = results.x(4,:);
    t = results.t_grid;
else
    q1 = results.extended.x(1,:);
    q2 = results.extended.x(2,:);
    dq1 = results.extended.x(3,:);
    dq2 = results.extended.x(4,:);
    t = results.extended.t_grid;
end

%% Prepare plotting
% Track
xtrack = 0:max(q1)/(length(q1)-1):max(q1);
h1 = plot_fc1(xtrack);
h2 = plot_fc2(xtrack);
% Positions
x1 = q1;
y1 = q2;

%% Plot results
% Cartesian
f1 = figure('Name', 'Position plot');
plot(x1,y1);
hold on
plot(xtrack,h1, LineStyle="--", Color='black')
plot(xtrack,h2, LineStyle="--", Color='black')
axis equal
grid on
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
legend('$q$', '$f_{c1}$', '$f_{c2}$','interpreter','latex', 'location', 'northeast')
axis([0 max(x1) min(y1) max(y1)+5])
% States
f2 = figure('Name', 'Generalized coordinates');
subplot(211)
plot(t,q1);
grid on
xlabel('$t$','interpreter','latex');
ylabel('$q_1$','interpreter','latex');
subplot(212)
plot(t,q2);
grid on
xlabel('$t$','interpreter','latex');
ylabel('$q_2$','interpreter','latex');
% Speeds
f3 = figure('Name', 'Generalized speeds');
subplot(211)
plot(t,dq1);
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v_1$','interpreter','latex');
subplot(212)
plot(t,dq2);
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v_2$','interpreter','latex');
% Speed of time
if ~b_ext
    f4 = figure('Name', 'Speed of time');
    subplot(121)
    plot(results.t_grid,t)
    hold on
    plot(results.t_grid,results.t_grid,'k--')
    grid on
    xlabel('$\tau$','interpreter','latex');
    ylabel('$t$','interpreter','latex');
    subplot(122)
    stairs(results.s_sot)
    grid on
    xlabel('simulation step','interpreter','latex');
    ylabel('$s$','interpreter','latex');
end
% Time lapse
f6 = figure('Name', 'Time lapse');
% plot(x1,y1, 'b');
hold on
plot(xtrack,h1, LineStyle="--", Color='black')
area(xtrack,h1, 'FaceColor', 'k', 'FaceAlpha', 0.2, 'EdgeAlpha', 0)
plot(xtrack,h2, LineStyle=":", Color='black')
area(xtrack,h2, 'FaceColor', 'k', 'FaceAlpha', 0.2, 'EdgeAlpha', 0)
px = [x1];
py = [y1];
axis([0 max(x1) min(y1) max(y1)+5])
for i = 1:5:size(px,2)
    hplot = scatter(px(:,i),py(:,i), 75, 'blue', 'filled');
    hold on
end
scatter(px(:,end),py(:,end), 75, 'blue', 'filled');
ylabel('$q_1$', 'interpreter','latex')
xlabel('$q_2$', 'interpreter','latex')
legend('$f_{c1}$', '', '$f_{c2}$', '', '$q(t)$', 'interpreter','latex', 'location', 'northeast')

%% Ski jump system
function [model, plot_fc1, plot_fc2] = ski_model_sim_nosnoc()
    import casadi.*
    model = NosnocModel();
    
    %% Model definition with slack
    % Physical parameters
    m1 = 80;
    g = 9.81;
    e = 1e-3;
    % States
    nq = 2;
    q = SX.sym('q',nq); 
    v = SX.sym('v',nq);
    model.x = [q;v];
    model.e = 0; % Coeff of restitution; 0 == inelastic
    p1 = -0.0001524;
    p2 = -0.002149;
    p3 = 0.4162;
    p4 = -9.154;
    p5 = 99.15;
    h1 = p1*q(1).^4+p2*q(1).^3+p3*q(1).^2+p4*q(1)+p5;
    plot_fc1 = @(q1) p1*q1.^4+p2*q1.^3+p3*q1.^2+p4*q1+p5; % function handle to plot contact surface
    h2 = -0.5*q(1)+45;
    plot_fc2 = @(q1) -0.5*q1+45; % function handle to plot contact surface
    model.x0 = [10; plot_fc1(10)+e; 0; 0];
    model.dims.n_dim_contact = 2;
    M = diag([m1 m1]);
    model.M = M;
    model.f_c = [q(2)-h1; q(2)-h2];
    model.f_v = model.M*[0; -g]; % Dynamics
    model.a_n = 100*g; % For auxiliary dynamical system
end