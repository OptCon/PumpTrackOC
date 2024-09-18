% Pump track benchmark system simulation

clear all;
clc;
close all

% REQUIRES CASADI AND NOSNOC IN MATLAB PATH
% Publication results obtained with: casadi-v3.5.5, nosnoc_v050
addpath('YOURPATH\casadi-windows-matlabR2016a-v3.5.5')
addpath('YOURPATH\nosnoc-main_v050\src')

import casadi.*

% Simulation or Optimization
mode = 'simulation';        % Input u == 0
% mode = 'optimization';    % Solving OCP 

%% Flags
b_ext = 0;

%% Parameters
N_FE = 3;
N_stages = 30;
T_sim = 2;
T0 = 1;
maxIter = 5e3;
q1_target = 2*pi;
R_U = 0.01;
P_T = 1000;
x0 = [0; 0; 0.4368; 22/3.6; 0; 0; 0]; % (smin+smax)/2

%% Initialize nosnoc
problem_options = NosnocProblemOptions();
solver_options = NosnocSolverOptions();

%% Settings
problem_options.time_freezing = 1;
problem_options.time_freezing_inelastic = 1;
problem_options.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
problem_options.N_finite_elements = N_FE;
solver_options.print_level = 3;  
solver_options.opts_casadi_nlp.ipopt.max_iter = maxIter;

if isequal(mode, 'simulation')  
    problem_options.stagewise_clock_constraint = 0;  
    problem_options.T_sim = T_sim;
    problem_options.N_sim = N_stages;
elseif isequal(mode, 'optimization')
    problem_options.time_optimal_problem = 1;
    problem_options.T = T0;
    problem_options.N_stages = N_stages;
else
    disp('Error: Mode unknown')
end

%% Initialize model
if isequal(mode, 'simulation')
    [model, plot_fc] = pump_model_sim_nosnoc(x0);
elseif isequal(mode, 'optimization')
    [model, plot_fc] = pump_model_opt_nosnoc(x0, q1_target, R_U, P_T);
else
    disp('Error: Mode unknown')
end

%% Solve
if isequal(mode, 'simulation')
    integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
    [results,stats] = integrator.solve();
elseif isequal(mode, 'optimization')
    % Solve OCP via nosnoc
    mpcc = NosnocMPCC(problem_options, model);
    solver = NosnocSolver(mpcc, solver_options);
    [results,stats] = solver.solve();
else
    disp('Error: Mode unknown') 
end

%% Unfold results
if ~b_ext
    q1 = results.x(1,:);
    q2 = results.x(2,:);
    q3 = results.x(3,:);
    dq1 = results.x(4,:);
    dq2 = results.x(5,:);
    dq3 = results.x(6,:);
    t = results.x(7,:);
else
    q1 = results.extended.x(1,:);
    q2 = results.extended.x(2,:);
    q3 = results.extended.x(3,:);
    dq1 = results.extended.x(4,:);
    dq2 = results.extended.x(5,:);
    dq3 = results.extended.x(6,:);
    t = results.extended.x(7,:);
end

%% Prepare plotting
% Track
h = plot_fc(q1);
% Cartesian coordinates
x1 = q1;
y1 = h+q2; % m_b
x2 = q1;
y2 = h+q2+q3; % m_r
% Horizontal speed
v1 = dq1;

%% Plot results

% Cartesian
f1 = figure('Name', 'Position in yz-plane');
plot(x1(1:2:end),y1(1:2:end), LineStyle="--", Color='blue');
hold on
plot(x2(1:2:end),y2(1:2:end), LineStyle="--", Color='blue');
area(x1(1:2:end),h(1:2:end), 'FaceColor', 'k', 'FaceAlpha', 0.2, 'EdgeAlpha', 0)
for i = 1:10:size(x1,2)
    scatter(x1(:,i),y1(:,i), 50, 'blue', 'filled');
    scatter(x2(:,i),y2(:,i), 50, 'blue', 'filled');
    line([x1(:,i),x2(:,i)], [y1(:,i),y2(:,i)], 'Color', 'blue');
    hold on
end
grid on
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
axis([0 max(q1) 0 1.25])

% States
f2 = figure('Name', 'Generalized coordinates');
subplot(311)
plot(t,q1);
axis([0 max(t) min(q1) max(q1)])
grid on
ylabel('$q_1$\,[m]','interpreter','latex');
subplot(312)
plot(t,q2);
axis([0 max(t) min(q2) max(q2)])
grid on
ylabel('$q_2$\,[m]','interpreter','latex');
subplot(313)
plot(t,q3);
yl = yline([0.59559 0.27803],'--',{'Max','Min'}, 'LabelHorizontalAlignment','left');
axis([0 max(t) 0.25 0.7])
grid on
xlabel('$t$\,[s]','interpreter','latex');
ylabel('$q_3$\,[m]','interpreter','latex');
% Speeds
f3 = figure('Name', 'Generalized speeds');
subplot(311)
plot(t,dq1);
axis([0 max(t) min(dq1) max(dq1)])
grid on
ylabel('$v_1\,[m/s]$','interpreter','latex');
subplot(312)
plot(t,dq2);
axis([0 max(t) min(dq2) max(dq2)])
grid on
ylabel('$v_2$\,[m/s]','interpreter','latex');
subplot(313)
plot(t,dq3);
grid on
xlabel('$t\,[s]$','interpreter','latex');
ylabel('$v_3$\,[m/s]','interpreter','latex');
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
% Input
if isequal(mode, 'optimization')
    u_opt = results.u;
    f5 = figure('Name', 'Input');
    stairs(results.t_grid_u,[u_opt,[nan]]);
%         yline([30.15 -8.66],'--',{'Max','Min'})
%         axis([0 max(results.t_grid_u) -10 33])
    grid on
    xlabel('$t$\,[s]','interpreter','latex');
    ylabel('$u$\,[N]','interpreter','latex');
end

%% Pump track benchmark system dynamics

function [model, plot_fc] = pump_model_sim_nosnoc(x0)
    import casadi.*
    model = NosnocModel();
    
    %% Model definition with slack
    % Physical parameters
    g = 9.81;
    m1 = 15; % mass bike
    m2 = 75; % mass rider
    % States
    nq = 3;
    q = SX.sym('q',nq); 
    v = SX.sym('v',nq);
    t = SX.sym('t',1);
    u = 0; % Input in simulation bugged
    model.x = [q;v;t];
    model.x0 = x0;
    model.e = 0; % Coeff of restitution; 0 == inelastic
    k_n1 = 50;
    h = -0.2*cos(2*q(1))+0.2;
    dh = 0.4*sin(2*q(1));
    ddh = 0.8*cos(2*q(1));
    plot_fc = @(q1) -0.2*cos(2*q1)+0.2; % function handle to plot contact surface
    model.c = [q(2); v(2)];
    model.S = [1 1; 1 -1; -1 1; -1 -1];
    f_ode = [   v(1);... % Unconstrained dynamics
                v(2);...
                v(3);...
                0;...
                - (m1*ddh*v(1)^2+u+g*m1)/m1;...
                (u*(m1+m2))/(m1*m2);...
                1]; % Time
    f_aux = k_n1.*[   0;... % Auxiliary dynamics
                0;...
                0;...
                -dh;...
                1;...
                0;...
                0]; % Time
    model.F = [f_ode f_ode f_ode f_aux]; % In matrix form
end


function [model, plot_fc] = pump_model_opt_nosnoc(x0, q1_target, R_U, P_T)
    import casadi.*
    model = NosnocModel();
    %% Model definition with slack
    % Physical parameters
    g = 9.81;
    m1 = 15; % mass bike
    m2 = 75; % mass rider
    % States
    nq = 3;
    q = SX.sym('q',nq); 
    v = SX.sym('v',nq);
    t = SX.sym('t',1);
    x = [q;v;t];
    model.x = x;
    model.ubx = [inf; inf; 0.59559; inf; inf; inf; inf]; % Results: at-paper
    model.lbx = [-inf; -inf; 0.27803; -inf; -inf; -inf; -inf]; % Results: at-paper
    model.x0 = x0;
    % Input
    u = SX.sym('u',1);
    model.u = u;
    model.lbu = -8.66*(m1*m2)/(m1+m2); % dds_min = -8.66
    model.ubu = 30.15*(m1*m2)/(m1+m2); % dds_max = 30.15
    model.u0 = 0;
    % Parameters
    model.e = 0; % Coeff of restitution; 0 == inelastic
    k_n = 50;
    h = -0.2*cos(2*q(1))+0.2;
    dh = 0.4*sin(2*q(1));
    ddh = 0.8*cos(2*q(1));
    plot_fc = @(q1) -0.2*cos(2*q1)+0.2; % function handle to plot contact surface
    % Switching functions and sign matrix
    model.c = [q(2); v(2)];
    model.S = [1 1; 1 -1; -1 1; -1 -1];
    % Dynamics
    f_ode = [       v(1);... % Unconstrained dynamics
                    v(2);...
                    v(3);...
                    0;...
                    -(m1*ddh*v(1)^2+u+g*m1)/m1;...
                    (u*(m1+m2))/(m1*m2);...
                    1]; % Time
    f_aux = k_n.*[  0;... % Auxiliary dynamics
                    0;...
                    0;...
                    -dh;...
                    1;...
                    0;...
                    0]; % Time
    model.F = [f_ode f_ode f_ode f_aux]; % In matrix form
    % Objective
%     model.f_q = R_U*(u'*u); % Objective term
%     model.f_q_T = P_T*(q(1)-q1_target).^2; % Terminal cost
    model.g_terminal = [q(1)-q1_target]; % Terminal constraint
end
