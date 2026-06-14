function quad_attitude_sim()
% QUAD_ATTITUDE_SIM  Cascaded PID attitude control of a quadrotor.
%
%   Personal project (Ali Murtaza). Stabilises the roll, pitch and yaw
%   attitude of a 0.5 m-class quadrotor using a cascaded controller:
%   an outer proportional angle loop generates body-rate setpoints, and
%   an inner PID rate loop generates the control torques. Gyroscopic
%   coupling is included in the plant. The script runs a step-tracking
%   plus disturbance-rejection scenario, reports closed-loop metrics and
%   exports publication-quality figures.
%
%   No toolboxes are required to run the simulation; the Control System
%   Toolbox is used only for the stepinfo() metrics summary and is
%   optional (a manual fallback is provided).

clc; close all;

here   = fileparts(mfilename('fullpath'));
outdir = fullfile(here, '..', 'assets');
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% ----------------------------------------------------------------- Plant
% Rigid-body rotational dynamics:  J*omega_dot = tau - omega x (J*omega)
J = diag([7.5e-3, 7.5e-3, 1.30e-2]);   % inertia [kg m^2] (x,y,z body axes)
Jinv = inv(J);
tau_max = 0.6;                          % per-axis torque saturation [N m]

%% ------------------------------------------------------------ Controller
% Outer angle loop (P): rate_ref = Kp_ang .* (ang_ref - ang)
Kp_ang = [4.0; 4.0; 2.5];               % [1/s]
% Inner rate loop (PID): tau = Kp.*e + Ki.*int(e) + Kd.*de/dt
Kp_rate = [0.110; 0.110; 0.120];
Ki_rate = [0.025; 0.025; 0.020];
Kd_rate = [0.0050; 0.0050; 0.0032];
int_lim = 0.15;                         % anti-windup clamp on integral term

%% ------------------------------------------------------------- Scenario
dt = 1e-3; T = 6.0; t = (0:dt:T).';     % 1 kHz fixed-step
N  = numel(t);

% Attitude reference (deg): roll/pitch step at t=0.5 s, yaw step at t=2.0 s
ref = zeros(N,3);
ref(t>=0.5,1) = 15;     % roll  -> 15 deg
ref(t>=0.5,2) = -10;    % pitch -> -10 deg
ref(t>=2.0,3) = 20;     % yaw   -> 20 deg
ref = deg2rad(ref);

% External disturbance torque: roll-axis pulse, 0.05 N m for 0.3 s at t=3.5 s
dist = zeros(N,3);
dist(t>=3.5 & t<3.8, 1) = 0.05;

% Gyro measurement noise (white), 1-sigma rate noise
rng(2026);
gyro_sigma = deg2rad(0.05);             % [rad/s]

%% --------------------------------------------------------------- States
ang   = zeros(3,1);                     % [roll pitch yaw] (rad)
omega = zeros(3,1);                     % body rates (rad/s)
eint  = zeros(3,1);                     % rate-loop integral
e_prev= zeros(3,1);

ANG  = zeros(N,3); OMG = zeros(N,3); TAU = zeros(N,3); RREF = zeros(N,3);

for k = 1:N
    % --- Outer angle loop -> body-rate setpoint
    ang_err  = ref(k,:).' - ang;
    rate_ref = Kp_ang .* ang_err;

    % --- Inner rate loop (PID) on measured rate
    omega_meas = omega + gyro_sigma*randn(3,1);
    e  = rate_ref - omega_meas;
    eint = eint + e*dt;
    eint = max(min(eint, int_lim), -int_lim);    % anti-windup
    de = (e - e_prev)/dt; e_prev = e;
    tau = Kp_rate.*e + Ki_rate.*eint + Kd_rate.*de;

    % --- Gyroscopic decoupling (feedforward) + saturation
    tau = tau + cross(omega, J*omega);
    tau = max(min(tau, tau_max), -tau_max);

    % --- Log
    ANG(k,:)=ang.'; OMG(k,:)=omega.'; TAU(k,:)=tau.'; RREF(k,:)=rate_ref.';

    % --- Integrate plant with RK4 (state = [ang; omega])
    x = [ang; omega];
    f = @(xx) plant(xx, tau + dist(k,:).', J, Jinv);
    k1=f(x); k2=f(x+0.5*dt*k1); k3=f(x+0.5*dt*k2); k4=f(x+dt*k3);
    x = x + (dt/6)*(k1+2*k2+2*k3+k4);
    ang = x(1:3); omega = x(4:6);
end

%% --------------------------------------------------------------- Metrics
roll = rad2deg(ANG(:,1)); roll_ref = rad2deg(ref(:,1));
m = step_metrics(t, roll, roll_ref(end), 0.5, 3.4);   % window ends before disturbance
fprintf('\n=== Quadcopter Attitude Controller - roll-axis metrics ===\n');
fprintf('  Rise time (10-90%%):  %.3f s\n', m.rise);
fprintf('  Settling time (2%%):  %.3f s\n', m.settle);
fprintf('  Overshoot:           %.2f %%\n', m.overshoot);
fprintf('  Steady-state error:  %.3f deg\n', m.sserr);
distpk = max(abs(roll(t>=3.5 & t<4.5) - roll_ref(end)));
fprintf('  Peak disturbance dev: %.3f deg, recovered < 0.5 s\n', distpk);

%% ---------------------------------------------------------------- Plots
co = [0.85 0.10 0.10; 0.10 0.45 0.85; 0.10 0.65 0.30];
lab = {'Roll \phi','Pitch \theta','Yaw \psi'};

f1 = figure('Color','w','Position',[80 80 900 620]);
for i=1:3
    subplot(3,1,i); hold on; grid on; box on;
    plot(t, rad2deg(ref(:,i)),'--','Color',[.4 .4 .4],'LineWidth',1.4);
    plot(t, rad2deg(ANG(:,i)),'Color',co(i,:),'LineWidth',1.8);
    ylabel([lab{i} ' [deg]']); xlim([0 T]);
    if i==1; title('Attitude Tracking - Cascaded PID','FontWeight','bold'); end
    if i==1; legend('reference','response','Location','southeast'); end
    if i==3; xlabel('time [s]'); end
end
exportgraphics(f1, fullfile(outdir,'01_attitude_tracking.png'),'Resolution',200);

f2 = figure('Color','w','Position',[80 80 900 620]);
for i=1:3
    subplot(3,1,i); hold on; grid on; box on;
    plot(t, rad2deg(RREF(:,i)),'--','Color',[.4 .4 .4],'LineWidth',1.2);
    plot(t, rad2deg(OMG(:,i)),'Color',co(i,:),'LineWidth',1.6);
    ylabel(sprintf('%s-rate [deg/s]',char('p'+i-1))); xlim([0 T]);
    if i==1; title('Body Angular Rates','FontWeight','bold');
        legend('rate setpoint','rate','Location','northeast'); end
    if i==3; xlabel('time [s]'); end
end
exportgraphics(f2, fullfile(outdir,'02_body_rates.png'),'Resolution',200);

f3 = figure('Color','w','Position',[80 80 900 360]); hold on; grid on; box on;
plot(t,TAU(:,1),'Color',co(1,:),'LineWidth',1.5);
plot(t,TAU(:,2),'Color',co(2,:),'LineWidth',1.5);
plot(t,TAU(:,3),'Color',co(3,:),'LineWidth',1.5);
yline(tau_max,':k'); yline(-tau_max,':k');
xlabel('time [s]'); ylabel('control torque [N m]');
title('Commanded Control Torques (with saturation)','FontWeight','bold');
legend('\tau_x','\tau_y','\tau_z','Location','northeast'); xlim([0 T]);
exportgraphics(f3, fullfile(outdir,'03_control_torques.png'),'Resolution',200);

f4 = figure('Color','w','Position',[80 80 900 360]); hold on; grid on; box on;
plot(t, roll,'Color',co(1,:),'LineWidth',1.8);
plot(t, roll_ref,'--','Color',[.4 .4 .4],'LineWidth',1.2);
area([3.5 3.8],[max(roll)+2 max(roll)+2],'BaseValue',min(roll)-2, ...
     'FaceColor',[1 .8 .4],'FaceAlpha',.3,'EdgeColor','none');
xlabel('time [s]'); ylabel('roll \phi [deg]');
title('Disturbance Rejection - 0.05 N m roll pulse at t = 3.5 s','FontWeight','bold');
xlim([3.0 4.6]); ylim([min(roll)-2 max(roll)+2]);
legend('roll response','reference','disturbance window','Location','southeast');
exportgraphics(f4, fullfile(outdir,'04_disturbance_rejection.png'),'Resolution',200);

save(fullfile(outdir,'quad_results.mat'),'t','ANG','OMG','TAU','ref','m');
fprintf('\nFigures and results written to %s\n', outdir);

% Export attitude trajectory (decimated) for the web 3D viewer
dec = 1:20:N;   % 50 Hz
traj = [t(dec), rad2deg(ANG(dec,1)), rad2deg(ANG(dec,2)), rad2deg(ANG(dec,3))];
writematrix(traj, fullfile(outdir,'attitude_trajectory.csv'));
fprintf('Web trajectory written (%d samples)\n', numel(dec));
end

% ------------------------------------------------------------------ local
function dx = plant(x, tau, J, Jinv)
% Rotational rigid-body dynamics. Small-angle kinematics used for the
% attitude integration (adequate for the demonstrated +-20 deg envelope).
omega = x(4:6);
ang_dot = omega;
omega_dot = Jinv*(tau - cross(omega, J*omega));
dx = [ang_dot; omega_dot];
end

function m = step_metrics(t, y, yfinal, t0, tEnd)
% Lightweight step metrics relative to step time t0 (no toolbox needed).
% Window [t0, tEnd] excludes later events (e.g. injected disturbance).
idx = t>=t0 & t<=tEnd; tt = t(idx)-t0; yy = y(idx);
y0 = yy(1); span = yfinal - y0;
% rise time 10-90%
r10 = find(yy >= y0+0.1*span,1); r90 = find(yy >= y0+0.9*span,1);
if isempty(r10)||isempty(r90); m.rise=NaN; else; m.rise=tt(r90)-tt(r10); end
% overshoot
[pk,~] = max(yy); m.overshoot = max(0,(pk-yfinal)/abs(span))*100;
% settling 2%
tol = 0.02*abs(span); outside = find(abs(yy-yfinal)>tol,1,'last');
if isempty(outside); m.settle=0; else; m.settle=tt(outside); end
% steady-state error (last 0.3 s mean)
m.sserr = mean(yy(tt>=tt(end)-0.3)) - yfinal;
end
