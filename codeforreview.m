clc;
clear all;
opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-14); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',100,'MaxIter',100,'TolFun',1e-10,'TolX',1e-14); % fsolve

% Initialize

% Constants
mu = 1;  % Gravitational parameter (AU^3/TU^2)
rf = 1.524;  % Mars final position in AU
r0 = 1.0;  % Earth initial position in AU

T = 70 / 58.13;  % Convert time to Time Units (TU)

% Initial Conditions
x0 = [r0, 0, 0, 1]; % Initial state: [r, theta, u, v]
xf = [rf, NaN, 0, sqrt(mu/rf) NaN 0 NaN NaN]; % Final state: [r, theta, u, v]

% Time
t0 = 0;
tf = T;

%% Numerical Solution (Min u)
lam0_guess = [0; 0; 0; 0];  % Initial guesses for costates
   
[lam0,fval] = fsolve(@cost_minU,lam0_guess,options,t0,tf,x0,xf,mu,opts_ode);
    
% Propagate Min-U Solution
[t_minU,X_minU] = ode45(@linosc,[t0 tf],[x0'; lam0],opts_ode,mu);
    time = t_minU(:,1);
    rr = X_minU(:,1);
    th = X_minU(:,2);
    uu = X_minU(:,3);
    vv = X_minU(:,4);
    lam_rr = X_minU(:,5);
    lam_th = X_minU(:,6);
    lam_uu = X_minU(:,7);
    lam_vv = X_minU(:,8);
    ur = -lam_uu;  % Radial control u_r
    utheta = -lam_vv;  % Tangential control u_theta




%% Functions
function Xdot = linosc(t,X,mu)

r   = X(1);       % Radial distance
theta = X(2);    % Angular position
u   = X(3);       % Radial control
v   = X(4);       % Tangential control
lam_r = X(5);     % Costate for r
lam_theta = X(6); % Costate for theta
lam_u = X(7);     % Costate for u
lam_v = X(8);     % Costate for v

% State Differential Equations
rdot = u;
thetadot = v / r;
udot = v^2 / r - mu / r^2 - lam_u;
vdot = -u * v / r - lam_v;

% Costate Differential Equations
lam_rdot = lam_theta * v / r^2 + lam_v * v^2 / r^2 - 2 * lam_u * mu / r^3 - lam_v * u * v / r;
lam_thetadot = 0;
lam_udot = -lam_r + lam_v * v / r
lam_vdot = -lam_theta / r - 2 * lam_u * v / r + lam_v * u / r

Xdot = [rdot; thetadot; udot; vdot; lam_rdot; lam_thetadot; lam_udot; lam_vdot];

end




function err = cost_minU(lam0_guess,t0,tf,x0,xf,mu,opts_ode)

[t,X] = ode45(@linosc,[t0 tf],[x0'; lam0_guess],opts_ode,mu);


err = [X(end,1) - xf(1), X(end,3) - 0, X(end,4) - xf(4), X(end,6) - 0];

end

% Performance Index
J = 0.5 * cumtrapz(time, ur.^2 + utheta.^2);

polarplot(th,rr);
title('Trajectory from Earth to Mars');

figure;

subplot(2, 1, 1);
plot(time, ur, 'r', 'LineWidth', 1.5); % Radial control u_r
title('Radial Control Input u_r');
xlabel('Time (TU)');
ylabel('u_r');

subplot(2, 1, 2);
plot(time, utheta, 'b', 'LineWidth', 1.5); % Tangential control u_theta
title('Tangential Control Input u_\theta');
xlabel('Time (TU)');
ylabel('u_\theta');

figure;
subplot(2,2,1);
plot(time, rr, 'b', 'LineWidth', 1.5);
title('Radial Distance (r)');
subplot(2,2,2);
plot(time, th, 'b', 'LineWidth', 1.5);
title('Angular Position (\theta)');
subplot(2,2,3);
plot(time, uu, 'b', 'LineWidth', 1.5);
title('Radial Velocity (u)');
subplot(2,2,4);
plot(time, vv, 'b', 'LineWidth', 1.5);
title('Tangential Velocity (v)');

figure;
subplot(2,2,1);
plot(time, lam_rr, 'b', 'LineWidth', 1.5);
title('Costate \lambda_r');
subplot(2,2,2);
plot(time, lam_th, 'b', 'LineWidth', 1.5);
title('Costate \lambda_\theta');
subplot(2,2,3);
plot(time, lam_uu, 'b', 'LineWidth', 1.5);
title('Costate \lambda_u');
subplot(2,2,4);
plot(time, lam_vv, 'b', 'LineWidth', 1.5);
title('Costate \lambda_v');

% Plot the cost function J versus time
figure;
plot(time, J, 'b', 'LineWidth', 1.5);
title('Cumulative Cost Function J(t)');
xlabel('Time (TU)');
ylabel('Cost J(t)');
J_total = J(end);  % This is the final value of the cumulative cost J(t)

% Plot the total cost vs. time of flight

T_values = [70, 100, 200, 360, 600];
J_values = zeros(size(T_values)); % To store the total cost for each time of flight
for i = 1:length(T_values)
    T = T_values(i) / 58.13;  % Convert time to Time Units (TU)
    
    % Initial Conditions
    x0 = [r0, 0, 0, 1]; % Initial state: [r, theta, u, v]
    xf = [rf, NaN, 0, sqrt(mu/rf), NaN, 0, NaN, NaN]; % Final state: [r, theta, u, v]
    
    % Time
    t0 = 0;
    tf = T;
    
    %% Numerical Solution (Min u)
    lam0_guess = [0; 0; 0; 0];  % Initial guesses for costates
       
    [lam0, fval] = fsolve(@cost_minU, lam0_guess, options, t0, tf, x0, xf, mu, opts_ode);
    
    % Propagate Min-U Solution
    [t_minU, X_minU] = ode45(@linosc, [t0 tf], [x0'; lam0], opts_ode, mu);
    time = t_minU(:, 1);
    rr = X_minU(:, 1);
    th = X_minU(:, 2);
    uu = X_minU(:, 3);
    vv = X_minU(:, 4);
    lam_rr = X_minU(:, 5);
    lam_th = X_minU(:, 6);
    lam_uu = X_minU(:, 7);
    lam_vv = X_minU(:, 8);
    ur = -lam_uu;  % Radial control u_r
    utheta = -lam_vv;  % Tangential control u_theta
    
    %% Performance Index
    J = 0.5 * cumtrapz(time, ur.^2 + utheta.^2);
    
    % Store the total cost for the current time of flight
    J_values(i) = J(end);
end

% Plot total cost J versus time of flight T
figure;
plot(T_values, J_values, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
title('Total Cost Function J versus Time of Flight');
xlabel('Time of Flight (T) [Time Units]');
ylabel('Total Cost J');
grid on;