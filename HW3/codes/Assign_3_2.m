%% Value approximation in LQR for Lecture 16 in 16-899 ACRL spring 2020
clc
clear all
close all
x0 = [0;0;0;0];
xg = [10;10;0;90*pi/180];
velocity = 0;
theta = 90 * pi / 180;

nx = 4;
nu = 2;

% linearized at xg
A =[1, 0,  0.1* cos(theta), 0.1 * -velocity * sin(theta); 
       0, 1, 0.1 * sin(theta), 0.1 * velocity * cos(theta); 
       0, 0,                      1,                                     0; 
       0, 0,                      0,                                    1];
B = [0, 0; 0, 0; 0.1, 0; 0, 0.1];
Q = eye(4);
R = 0.1 * eye(2);
S = 10 * eye(4);
N = 100;

P = cell(N);
P{N} = S;
W = 0.0001* eye(1);
w = normrnd(0,sqrt(W), [1,N]);

% Solving Riccati Equation
for i = N-1:-1:1
    P{i} = Q + A' * P{i+1} * A - A' * P{i+1} * B * inv(R + B' * P{i+1} * B) * B' * P{i+1} * A;    p11(i) = P{i}(1,1);
end


K = - inv(B' * P{1} * B + R) * B' * P{1} * A;
p.A = A; p.B = B; p.Q = Q; p.R = R;
p.terminal = 0.01; p.nx = nx; p.nu = nu; p.x0 = x0;
x = x0;
k = 0; 
k_max = 100;
traj_opt = zeros(1,nx*k_max);
cost = 0;

traj = zeros(4,1,k_max);
traj(:,:,1) = x0;
while ~terminal( x- xg, p) && k < k_max
    u = K*(x-xg);
%     [xnew, l] = dynamics(x, u, p, xg);

    l = ((x - xg)'*p.Q*(x-xg) + u'*p.R*u)/2;
    xnew = x + [x(3)*cos(x(4)); x(3)*sin(x(4)); u(1); u(2)]*0.1 + w(i);
    
    cost = cost + l;
    x = xnew;
    k = k+1;
    traj_opt(1,nx*(k-1)+1:nx*k) = x;
    traj(:,:,k) = x;
end

W_gt = [Q+A'*P{1}*A  A'*P{1}*B; 
             B'*P{1}*A       B'*P{1}*B+R];

figure(1)
plot(squeeze(traj(1,1,:)))
hold on 
plot(squeeze(traj(2,1,:)))
hold on 
plot(squeeze(traj(3,1,:)))
hold on 
plot(squeeze(traj(4,1,:)))


%% Learning paramters
W_init = [1 1 2 1 2 3; 
               1 2 1 3 2 1;
               2 1 2 3 1 4;
               1 3 3 4 1 2;
               2 2 1 1 1 2;
               3 1 4 2 2 1];
n_ep = 100;
epsilon = 0.5;
delta = 0.1;
alpha = 0.01;
k_max = 20;
%% Monte Carlo
clc
close all
T = 0.1;
W = W_init; % Qudratic parameterization
% traj_mc = zeros(n_ep, nx*k_max);
rms_mc = zeros(1,n_ep);
figure(1);
clf;
trajectory = zeros(4,1,k_max);
for episode = 1:n_ep
    x = x0; 
    k = 1; 
    u = greedy(x, W, epsilon, xg);
    dW = zeros(nx+nu);
    x_list = zeros(4,1,k_max);
    x_list(:,:,1) = x;
    u_list = zeros(2,1,k_max);
    u_list(:,:,1) = u;
    l_list = [];
    while ~terminal(x-xg, p) && k <= k_max
%         [x, l] = dynamics(x, u, p, xg);
        l = ((x - xg)'*p.Q*(x-xg) + u'*p.R*u)/2; %run-time cost
        x = x + [x(3)*cos(x(4)); x(3)*sin(x(4)); u(1); u(2)]*0.1 ;
        u = greedy(x, W, epsilon,xg);
%         x_list = [x_list;x];
%         u_list = [u_list;u];
        x_list(:,:,k+2) = x;
        u_list(:,:,k+2) = u;
        l_list = [l_list;l];
        k = k+1;
%         traj_mc(episode, nx*(k-1)+1:nx*k) = x;
        trajectory(:,:,k) = x;
    end
    subplot(411)
    plot(squeeze(trajectory(1,1,:)))
    hold on
    subplot(412)
    plot(squeeze(trajectory(2,1,:)))
    hold on
    subplot(413)
    plot(squeeze(trajectory(3,1,:)))
    hold on
    subplot(414)
    plot(squeeze(trajectory(4,1,:)))
    hold on
    
    for kk = 1: k
        G = sum(l_list(kk:end));
        x = x_list(:,:,kk); 
        u = u_list(:,:,kk);
        f = [x(1) + x(3) *cos(x(4))*T - xg(1);
               x(2) + x(3) *sin(x(4))*T - xg(2);
               x(3) - xg(3);
               x(4) - xg(4)
             ];
        dW = dW - alpha * (G - [f;u]'*W*[f;u]/2) * 0.5*[f;u]*[f;u]';
    end
    W = W + dW;
    rms_mc(episode) = norm(W - W_gt);
end
figure(2)
plot(rms_mc(1:9))
%% SARSA
close all
W_init = [1 1 2 1 2 1;
               1 3 4 1 2 2;
               2 4 5 6 1 3;
               1 1 6 2 3 3;
               2 2 1 3 2 1;
               1 2 3 3 1 0 ];
alpha = 0.0001;
W = W_init; % Qudratic parameterization
% traj_sarsa = zeros(n_ep,nx*k_max);
rms_sarsa = zeros(1,n_ep);
% subplot(312);hold on;
figure(1);
clf;
trajectory = zeros(4,1,k_max);
for episode = 1:n_ep
    x = x0; k = 0; 
    u = greedy(x, W, epsilon, xg);
    while ~terminal(x-xg, p) && k < k_max
%         [xnew, l] = dynamics(x, u, p, xg);
        l = -((x - xg)'*p.Q*(x-xg) + u'*p.R*u)/2;
        xnew = x + [x(3)*cos(x(4)); x(3)*sin(x(4)); u(1); u(2)]*0.1 + w(i);
       T = 0.1;
        unew = greedy(xnew, W, epsilon, xg);
        f = [x(1) + x(3) *cos(x(4))*T - xg(1);
               x(2) + x(3) *sin(x(4))*T - xg(2);
               x(3) - xg(3);
               x(4) - xg(4)
             ];
        W = W + alpha * (l + delta * [f;unew]'*W*[f;unew]/2 - [f;u]'*W*[f;u]/2) * 0.5*[f;u]*[f;u]';
        x = xnew;
        u = unew;
        k = k+1;
%         traj_sarsa(episode,nx*(k-1)+1:nx*k) = x;
        trajectory(:,:,i) = x;
    end
    subplot(411)
    plot(squeeze(trajectory(1,1,:)))
    hold on
    subplot(412)
    plot(squeeze(trajectory(2,1,:)))
    hold on
    subplot(413)
    plot(squeeze(trajectory(3,1,:)))
    hold on
    subplot(414)
    plot(squeeze(trajectory(4,1,:)))
    hold on
%     plot(0:k,[x0 traj_sarsa(episode,1:nx*k)],'color',[0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep])
    rms_sarsa(episode) = norm(W - W_gt);
end
% plot(0:10,[x0 traj_opt(1,1:nx*10)],'--k')
% title("Sarsa")
% box on;
%% Q-Learning
W_init = [1600 1000  2000 1500 2500 2000; 
               1000 1800 1000 1300 2100 1200;
               2000 1000 1500 2000 3002 1000;
               1500 1300 2000 1400 1100 1000;
               2500 2100 3200 1001 1400 1200;
               2000 1200 1200 1000 1200 1500];
W = W_init; % Qudratic parameterization
% traj_q = zeros(n_ep,nx*k_max);
rms_q = zeros(1,n_ep);
% subplot(313);hold on;
figure(1);
clf;
trajectory = zeros(4,1,k_max);
alpha = 0.0001;
for episode = 1:n_ep
    x = x0; k = 0;
    while ~terminal(xg-x, p) && k < k_max
        u = greedy(x, W, epsilon,xg);
%         [xnew, l] = dynamics(x, u, p);
        l = ((xg - x)'*p.Q*(xg-x) + u'*p.R*u)/2;
        xnew = x + [x(3)*cos(x(4)); x(3)*sin(x(4)); u(1); u(2)]*0.1 + w(i);
  
        W = W + alpha * (l + delta * (xg - xnew)'*min_Q(W,p)*(xg - xnew)/2 - [xg-x;u]'*W*[xg-x;u]/2) * 0.5*[xg-x;u]*[xg-x;u]'
        x = xnew;
        k = k+1;
%         traj_q(episode,nx*(k-1)+1:nx*k) = x;
        trajectory(:,:,i) = x;
    end
    subplot(411)
    plot(squeeze(trajectory(1,1,:)))
    hold on
    subplot(412)
    plot(squeeze(trajectory(2,1,:)))
    hold on
    subplot(413)
    plot(squeeze(trajectory(3,1,:)))
    hold on
    subplot(414)
    plot(squeeze(trajectory(4,1,:)))
    hold on
%     plot(0:k,[x0 traj_q(episode,1:nx*k)],'color',[0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep])
    rms_q(episode) = norm(W - W_gt);
end
%% Functions
function V = min_Q(W,p)
nx = p.nx;
W11= W(1:nx,1:nx);
W12= W(1:nx,nx+1:end);
W21= W(nx+1:end, 1:nx);
W22= W(nx+1:end, nx+1:end);
V = W11 - W12 * inv(W22) * W21;
end
function [xnew, cost] = dynamics(x,u,p,xg)
cost = ((x - xg)'*p.Q*(x-xg) + u'*p.R*u)/2;
xnew = p.A*x + p.B*u + [xg(4)*xg(3)*sin(xg(4))*0.1; - xg(4)*xg(3)*cos(xg(4))*0.1;0;0 ];
end
% 
function flag = terminal(x,p)
    if norm(x) < p.terminal
        flag = true;
    else
        flag = false;
    end
end

function u = greedy(x, W, epsilon, xg)
    T = 0.1;
    sample = rand(1);
    nx = length(x);
    nu = size(W,1)-nx;
    f = [x(1) + x(3) *cos(x(4))*T - xg(1);
           x(2) + x(3) *sin(x(4))*T - xg(2);
           x(3) - xg(3);
           x(4) - xg(4)
        ];
    u = - inv(W(nx+1:end, nx+1:end))*(W(nx+1:end, 1:nx) + W(1:nx, nx+1:end).')/2*f;
    if sample < epsilon
        u = -rand(nu,nx)*f;
    end
end