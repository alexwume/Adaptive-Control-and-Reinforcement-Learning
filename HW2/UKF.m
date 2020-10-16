%% problem 4.3
close all
clear all
%number of steps
N = 1000;

Q = eye(2);
R = 0.1;


A = cell(N);
B = cell(N);

%parameters
a = 0.12;
b = 0.008;
c = 0.11;
a_hat = zeros(N,1);
b_hat = zeros(N,1);
c_hat = zeros(N,1);
a_hat(1,1) = 0.1;
b_hat(1,1) = 0.005;
c_hat(1,1) = 0.1;

%states
x = zeros(2,1,N);
x(:,:,1) = [2;0];
x_hat = zeros(2,1,N);
x_hat(:,:,1) = [3;0];
x_extended = zeros(5,1,N);
x_extended(:,:,1) = [x_hat(:,:,1); a_hat(1); b_hat(1); c_hat(1)];

%outputs
y = zeros(2,1,N);
y(:,:,1) = x(:,:,1);

%actual system
A_a = [1 a; 0 1];
B_a = [b; c];

%noise
W = 0.0001* eye(1);
W_bar = diag([W, 0, 0, 0]);
V = 0.00001*eye(2);
w = normrnd(0,sqrt(W), [1,N]);
w_theta = [0.000;0.000;0.000];

itr = 30;
control = zeros(N,1);

Zk = blkdiag(0.001 * eye(2), w_theta(1), w_theta(2), w_theta(3));

for i = 2:N
    A{i -1} = [1 a_hat(i - 1); 0 1];
    B{i - 1}= [b_hat(i - 1) ; c_hat(i - 1)];

    P = [2 5; 1 10];
    
    % Solving Riccati Equation for LQR control
    for j = itr:-1:1
        P = Q + A{i - 1}' * P * A{i - 1} - A{i - 1} * P * B{i - 1} * inv(R + B{i - 1}' * P * B{i - 1}) * B{i - 1}' * P * A{i - 1};
    end
    control_gain = -inv(R + B{i - 1}' * P * B{i - 1}) * B{i - 1}' * P * A{i - 1};

    control(i - 1) = control_gain * x_hat(:,:,i - 1);


    B_bar = blkdiag(B{i - 1}, 1,1,1);
    C_bar = [1 0 0 0 0; 0 1 0 0 0];
         
    %UKF sampling
    sample_states = zeros(5,itr);
    sample_states(:,:) =mvnrnd(x_extended(:,:, i - 1), Zk,itr).'; %samples, size = 5 x 30
    
    sample_states_prop= zeros(5,itr);
    sample_states_prop(:,:) = [A{i - 1} * sample_states(1:2,:); sample_states(3:5,:)] + [B{i - 1} * control(i - 1); 0; 0; 0 ];
       
    %UKF Dynamic updates
    x_extended_prio = mean(sample_states_prop,2); 
    Mk = cov(sample_states_prop.') + B_bar *W_bar* B_bar.';

    
    %UKF Measurement updates
    sample_states_measurement = zeros(5,itr);
    sample_states_measurement(:,:) =mvnrnd(x_extended_prio, Mk, itr).'; %samples, size = 5 x 30
    
    tmp_y = C_bar * sample_states_measurement(:,:);
    mean_y = mean(tmp_y, 2);
%     sig_y_k = diag(var(tmp_y,0,2)) + V;
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
    
    sig_xy_k = [cov11(1,2) cov12(1,2); cov21(1,2) cov22(1,2); cov31(1,2) cov32(1,2); cov41(1,2) cov42(1,2); cov51(1,2) cov52(1,2)];
    
    %Parameters updates
    KF_gain = sig_xy_k * inv(sig_y_k);
    measurement_noise = (normrnd(0,sqrt(V))) * [1;1];
    
    %actual sys
    x(:,:,i) = A_a * x(:,:,i -1) + B_a * (control(i - 1));
    y(:,:,i) = x(:,:,i) + measurement_noise;
    
    x_extended(:,:,i) = x_extended_prio + KF_gain * (y(:,:,i) - mean_y);

    x_hat(:,:, i) = x_extended(1:2,:,i);
    a_hat(i) = x_extended(3,:,i);
    b_hat(i) = x_extended(4,:,i);
    c_hat(i) = x_extended(5,:,i);
    
    Zk = Mk - KF_gain * (sig_xy_k.');
    for k = 1:1:5
        for j = 1:1:5
                Zk(k,j) = abs(Zk(k,j));
        end
    end
    
end

figure(1);
plot(squeeze(x_hat(1,1,:)))
hold on
plot(squeeze(y(1,1,:)))

figure(2);
plot(squeeze(x_hat(2,1,:)))
hold on
plot(squeeze(y(2,1,:)))

figure(3);
plot(squeeze(a_hat(:)));
hold on;
plot(squeeze(b_hat(:)));
hold on;
plot(squeeze(c_hat(:)));
