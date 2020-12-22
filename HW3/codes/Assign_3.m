%% RLS
close all
clear all
clc
N = 100;
T = 0.1;

a = 0.1;
b = 0.1;

a_hat = zeros(N,1);
b_hat = zeros(N,1);

a_hat(1,1) = 0.5;
b_hat(1,1) = 0.5;


x = zeros(2,1,N);
x0 = [0;0];
xg = [10;10];

x(:,:,1) = x0;
x_hat = zeros(2,1,N);
x_hat(:,:,1) = x0;

xg_hat = zeros(2,1,N);
xg_hat(:,:,1) = [6;8];

G = zeros(2,2);
H = zeros(2,2);
lambda = 0.5;
error = zeros(2,1);

W = 0.0001* eye(1);
w = normrnd(0,sqrt(W), [1,N]);

for i = 2:N
    %actual system
    x(:,:,i) = dynamics(xg, x(:,:,i - 1), a, b, T) + [T; T] * w(i);
    
    %parameter updates
    x_hat(:,:,i) = dynamics(xg_hat(:,:,i - 1), x(:,:,i - 1), a, b, T) + [T; T] * w(i);
    error = x(:,:,i) - x_hat(:,:,i);
    G = [a*T b* T; a*T b*T];
%     H = G.'*G + lambda * H;
    xg_hat(:,:, i) = xg_hat(:,:,i - 1) + 500*G.'*error;


end

figure(1);
plot(squeeze(xg_hat(1,1,:)));
hold on;
plot(squeeze(xg_hat(2,1,:)));
legend("Estimated xg1","Estimated xg2")
figure(2);
plot(squeeze(x(1,1,:)));
hold on;
plot(squeeze(x(2,1,:)));
legend("x1","x2")


%% problem 1.2
close all
clear all
clc
N = 100;
T = 0.1;

a = 0.1;
b = 0.1;

a_hat = zeros(N,1);
b_hat = zeros(N,1);

a_hat(1) = 0.4;
b_hat(1) = 0.01;


x = zeros(2,1,N);
x0 = [0;0];
xg = [10;10];

x(:,:,1) = x0;
x_hat = zeros(2,1,N);
x_hat(:,:,1) = x0;

theta = zeros(2,1,N);
theta(:,:,1) = [a_hat(1), b_hat(1)];

G = zeros(2,2,N);
H = [85000 20200; 200000 500000];
lambda = 0.8;
error = zeros(2,1,N);

W = 0.0001* eye(1);
w = normrnd(0,sqrt(W), [1,N]);

for i = 2:N
    %actual system
    x(:,:,i) = dynamics(xg, x(:,:,i - 1), a, b, T) + [T; T] * w(i);
    
    %parameter updates
    x_hat(:,:,i) = dynamics(xg, x(:,:,i - 1), a_hat(i - 1), b_hat(i - 1), T) + [T; T] * w(i);
    error(:,:,i) = x(:,:,i) - x_hat(:,:,i);
    G(:,:,i - 1) = [(xg(1) - x(1,1,i - 1))*T (xg(2) - x(2,1,i - 1))*T;
                        (xg(1) - x(1,1,i - 1))*T (xg(2) - x(2,1,i - 1))*T];
    H = G(:,:,i - 1).'*G(:,:,i - 1) + lambda * H;
    theta(:,:, i) = theta(:,:,i - 1) + inv(H)*G(:,:,i - 1).'*error(:,:,i);
    a_hat(i) = theta(1,1,i);
    b_hat(i) = theta(2,1,i);


end

figure(1);
plot(a_hat);
hold on;
plot(b_hat);
legend("Estimated a","Estimated b")

figure(2);
plot(squeeze(x(1,1,:)));
hold on;
plot(squeeze(x(2,1,:)));
hold on
plot(squeeze(x_hat(1,1,:)));
hold on;
plot(squeeze(x_hat(2,1,:)));

figure(3)
plot(squeeze(error(1,1,:)));
hold on;
plot(squeeze(error(2,1,:)));

%% problem 1.3
close all
clear all
clc
N = 100;
T = 0.1;

a = 0.1;
b = 0.1;

a_hat = zeros(N,1);
b_hat = zeros(N,1);

a_hat(1) = 0.5;
b_hat(1) = 0.2;


x = zeros(2,1,N);
x0 = [0;0];
xg = [10;10];
xg_hat = zeros(2,1,N);
xg_hat(:,:,1) = [2;4];


x(:,:,1) = x0;
x_hat = zeros(2,1,N);
x_hat(:,:,1) = x0;

theta = zeros(4,1,N);
theta(:,:,1) = [a_hat(1), b_hat(1), xg_hat(1,1,1), xg_hat(2,1,1)];

G = zeros(2,4,N);
H = [30000 20000 60000 20000; 
        60000 80000 70000 70000; 
        10000 50300 60000 20000; 
        80000 10000 40000 30000];
lambda = 0.80;
error = zeros(2,1,N);

W = 0.0001* eye(1);
w = normrnd(0,sqrt(W), [1,N]);

for i = 2:N
    %actual system
    x(:,:,i) = dynamics(xg, x(:,:,i - 1), a, b, T) + [T; T] * w(i);
    
    %parameter updates
    x_hat(:,:,i) = dynamics(xg_hat(:,:,i - 1), x(:,:,i - 1), a_hat(i - 1), b_hat(i - 1), T) + [T; T] * w(i);
    error(:,:,i) = x(:,:,i) - x_hat(:,:,i);
    G(:,:,i - 1) = [(xg_hat(1,1,i - 1) - x(1,1,i - 1))*T (xg_hat(2,1,i - 1) - x(2,1,i - 1))*T a_hat(i - 1)*T b_hat(i - 1)* T; 
                        (xg_hat(1,1,i - 1) - x(1,1,i - 1))*T (xg_hat(2,1,i - 1) - x(2,1,i - 1))*T a_hat(i - 1)*T b_hat(i - 1)* T];
                    
    H = G(:,:,i - 1).'*G(:,:,i - 1) + lambda * H;
    
    theta(:,:, i) = theta(:,:,i - 1) + 0.5*pinv(H)*G(:,:,i - 1).'*error(:,:,i);
%     inv(H)*G(:,:,i - 1).'*error(:,:,i)
%     theta(:,:,i)
%     error(:,:,i)
    a_hat(i) = theta(1,1,i);
    b_hat(i) = theta(2,1,i);
    xg_hat(1,1,i) = theta(3,1,i);
    xg_hat(2,1,i) = theta(4,1,i);

end

figure(1);
plot(a_hat);
hold on;
plot(b_hat);
legend("Estimated a","Estimated b")

figure(2);
plot(squeeze(x(1,1,:)));
hold on;
plot(squeeze(x(2,1,:)));
hold on

figure(3)
plot(squeeze(xg(1,1,:)));
hold on
plot(squeeze(xg(2,1,:)));

%% EKF
close all
clear all
N = 100;
Q = eye(2);
R = 0.1;
A = cell(N);
B = cell(N);

a = 0.1;
b = 0.1;

xg = [10; 10];
a_hat = zeros(N,1);
b_hat = zeros(N,1);
xg_hat = zeros(2,1,N);
a_hat(1,1) = 1;
b_hat(1,1) = 1;
xg_hat(:,:,1) = [2;0];

x = zeros(2,1,N);
x(:,:,1) = [0;0];

x_extended = zeros(6,1,N);
x_extended(:,:,1) = [x(:,:,1); a_hat(1); b_hat(1);xg_hat(1,1,1); xg_hat(2,1,1)];

W = 0.0001* eye(1);
w = normrnd(0,sqrt(W), [1,N]);
a_w = normrnd(0,0, [1,N]);
b_w = normrnd(0,0, [1,N]);
a_w = zeros(N);
b_w = zeros(N);
xg1_w = zeros(N);
xg2_w = zeros(N);
% xg1_w = normrnd(0,0, [1,N]);
% xg2_w = normrnd(0,0, [1,N]);
V = 0.0000001*eye(2);

control = zeros(N,1);
Zk = 1 * eye(6);

actual_control = zeros(N, 1);
%actual system
A = [1 0; 0 1];
B = [0.1; 0.1];
T = 0.1;

for i = 2:N
    control(i - 1) = [a_hat(i - 1) b_hat(i - 1)] * (xg_hat(:,:,i - 1) - x(:,:,i - 1));
  
    xg1 = xg_hat(1,:,i - 1);
    xg2 = xg_hat(2,:,i - 1);
    x1 = x(1,:,i - 1);
    x2 = x(2,:,i - 1);
%     x1 = 10;
%     x2 = 10;
    A_bar =[
                    1 - T*a_hat(i - 1),  -T *b_hat(i - 1), T*xg1 - T*x1, T*xg2 - T*x2, T*a_hat(i - 1), T*b_hat(i-1); 
                    -T*a_hat(i - 1),  1 - T*b_hat(i - 1), T*xg1 - T*x1, T*xg2 - T*x2, T*a_hat(i - 1), T*b_hat(i-1); 
                    0, 0, 1, 0, 0, 0;
                    0, 0, 0, 1, 0, 0;
                    0, 0, 0, 0, 1, 0;
                    0, 0, 0, 0, 0, 1;
                ];
    B_bar = blkdiag(B, eye(4));
    C_bar = [1 0 0 0 0 0; 0 1 0 0 0 0];
    
%     x_extended_prio = [A * x(:,:,i - 1) ; x_extended(3:6,1,i - 1)] + [B * control(i - 1); 0; 0; 0;0 ] + B_bar * [w(i - 1); a_w(i - 1); b_w(i - 1); xg1_w(i - 1); xg2_w(i - 1)];
    x_extended_prio = A_bar * x_extended(:,:,i - 1) +  B_bar * [w(i - 1); a_w(i - 1); b_w(i - 1); xg1_w(i - 1); xg2_w(i - 1)];
    Mk = A_bar * Zk * A_bar.' + B_bar * diag([w(i - 1), a_w(i - 1), b_w(i - 1), xg1_w(i - 1), xg2_w(i - 1)]) * B_bar.';
    
    KF_gain = Mk * C_bar.' * pinv(C_bar * Mk * C_bar.' + V);

    %actual sys
    actual_control(i - 1) = [a b] * (xg - x(:,:,i - 1));
    x(:,:,i) = A * x(:,:,i -1) + B * actual_control(i - 1) + B*w(i - 1);
    y(:,:,i) = x(:,:,i);
    
    x_extended(:,:,i) = x_extended_prio + KF_gain * (y(:,:,i) - C_bar * x_extended_prio);
    a_hat(i) = x_extended(3,:,i);
    b_hat(i) = x_extended(4,:,i);
    xg_hat(:,:,i) = x_extended(5:6, : , i);
    x_extended(1:2,1,i) = x(:,:,i); 
    Zk = Mk - KF_gain * C_bar * Mk;
    
    
end

figure(1);
plot(squeeze(x(1,1,:)))
legend('x1');

figure(2);
plot(squeeze(xg_hat(1,1,:)))
hold on
plot(squeeze(xg_hat(2,1,:)))
hold on
legend('xg1','xg2');

figure(3);
plot(squeeze(a_hat(:)));
hold on;
plot(squeeze(b_hat(:)));
legend('a = 0.1', 'b = 0.1');


%% problem 4.3
close all
clear all
%number of steps
N = 100;

%parameters
a = 0.1;
b = 0.1;

a_hat = zeros(N,1);
b_hat = zeros(N,1);

a_hat(1,1) = 0.12;
b_hat(1,1) = 0.005;


%states
x = zeros(2,1,N);
x(:,:,1) = [2;0];
xg_hat = zeros(2,1,N);
xg_hat(:,:,1) = [2;0];
x_extended = zeros(6,1,N);
x_extended(:,:,1) = [x(:,:,1); a_hat(1); b_hat(1); xg_hat(:,:,1)];

%outputs
y = zeros(2,1,N);
y(:,:,1) = x(:,:,1);

T = 0.1;
%actual system
A = [1 0; 0 1];
B = [T; T];

%noise
W = 0.0001* eye(1);
W_bar = diag([W, 0, 0, 0, 0]);
V = 0.00001*eye(2);
w = normrnd(0,sqrt(W), [1,N]);
w_theta = [1.000;1.000;1.000;1.000];

itr = 30;
control = zeros(N,1);

Zk = blkdiag(0.001 * eye(2), w_theta(1), w_theta(2), w_theta(3),w_theta(4));
% Zk = 

for i = 2:N

    %estimation
    control(i - 1) = [a_hat(i - 1) b_hat(i - 1)] * (xg_hat(:,:,i - 1) - x(:,:,i - 1));
  
    xg1 = xg_hat(1,:,i - 1);
    xg2 = xg_hat(2,:,i - 1);
    x1 = x(1,:,i - 1);
    x2 = x(2,:,i - 1);

    B_bar = blkdiag(B, eye(4));
    C_bar = [1 0 0 0 0 0; 0 1 0 0 0 0];

    %UKF sampling
    sample_states = zeros(6,itr);
    sample_states(:,:) =mvnrnd(x_extended(:,:, i - 1), Zk,itr).'; %samples, size = 6 x 30
    
    sample_states_prop= zeros(6,itr);
    sample_states_prop(:,:) = [A * sample_states(1:2,:); sample_states(3:6,:)] + [B * control(i - 1); 0; 0; 0; 0];
       
    %UKF Dynamic updates
    x_extended_prio = mean(sample_states_prop,2); 
    Mk = cov(sample_states_prop.') + B_bar *W_bar* B_bar.';

    
    %UKF Measurement updates
    sample_states_measurement = zeros(6,itr);
    sample_states_measurement(:,:) =mvnrnd(x_extended_prio, Mk, itr).'; %samples, size = 5 x 30
    
    tmp_y = C_bar * sample_states_measurement(:,:);
    mean_y = mean(tmp_y, 2);
    sig_y_k = cov(tmp_y.') + V;

    %constructing the covariance matrix of sample_states_prop and tmp_y
    cov11 = cov(sample_states_measurement(1,:),tmp_y(1,:));
    cov12 = cov(sample_states_measurement(1,:),tmp_y(2,:));
    cov21 = cov(sample_states_measurement(2,:),tmp_y(1,:));
    cov22 = cov(sample_states_measurement(2,:),tmp_y(2,:));
    cov31 = cov(sample_states_measurement(3,:),tmp_y(1,:));
    cov32 = cov(sample_states_measurement(3,:),tmp_y(2,:));
    cov41 = cov(sample_states_measurement(4,:),tmp_y(1,:));
    cov42 = cov(sample_states_measurement(4,:),tmp_y(2,:));
    cov51 = cov(sample_states_measurement(5,:),tmp_y(1,:));
    cov52 = cov(sample_states_measurement(5,:),tmp_y(2,:));
    cov61 = cov(sample_states_measurement(6,:),tmp_y(1,:));
    cov62 = cov(sample_states_measurement(6,:),tmp_y(2,:));
    
    
    sig_xy_k = [cov11(1,2) cov12(1,2); cov21(1,2) cov22(1,2); cov31(1,2) cov32(1,2); cov41(1,2) cov42(1,2); cov51(1,2) cov52(1,2);cov61(1,2) cov62(1,2)];
    
    %Parameters updates
    KF_gain = sig_xy_k * pinv(sig_y_k);
    measurement_noise = (normrnd(0,0)) * [1;1];
    
    %actual sys
    x(:,:,i) = A * x(:,:,i -1) + B * (control(i - 1));
    y(:,:,i) = x(:,:,i);
    
    x_extended(:,:,i) = x_extended_prio + KF_gain * (y(:,:,i) - mean_y);

    x_hat(:,:, i) = x_extended(1:2,:,i);
    a_hat(i) = x_extended(3,:,i);
    b_hat(i) = x_extended(4,:,i);
    xg_hat(:,:,i) = x_extended(5:6,:,i);
    
    Zk = Mk - KF_gain * (sig_xy_k.');
%     for k = 1:1:6
%         for j = 1:1:6
%                 Zk(k,j) = abs(Zk(k,j));
%         end
%     end
    
end

% figure(1);
% plot(squeeze(x_hat(1,1,:)))
% hold on
% plot(squeeze(y(1,1,:)))
% 
% figure(2);
% plot(squeeze(x_hat(2,1,:)))
% hold on
% plot(squeeze(y(2,1,:)))
% 
% figure(3);
% plot(squeeze(a_hat(:)));
% hold on;
% plot(squeeze(b_hat(:)));
% hold on;
% plot(squeeze(c_hat(:)));
% x_extended = zeros(5,1,N);
% x_extended(:,:,1) = [x(:,:,1); a_hat(1); b_hat(1); xg1_hat(1);xg2_hat(1)];

%% functions
function f = dynamics(xg, xk, a, b, T)
     f = zeros(2,1);
     u1 = a * xg(1) * T - a * xk(1) * T + b * xg(2) * T - b * xg(2) * T;
%      u1 = control_constraint(u1, 5, -5);
     u2 = a * xg(1) * T - a * xk(1) * T + b * xg(2) * T - b * xg(2) * T;
%      u2 = control_constraint(u2, 5, -5);
     f(1)= xk(1) + u1;
     f(2)= xk(2) + u2;
end

function u = control_constraint(input, upper, lower)
     if input > upper
         u = upper;
     elseif input < lower
         u = lower;
     else
         u = input;
     end
end


