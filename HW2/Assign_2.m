%% problem 1.2
x_bar = 5;
u_bar = 1;
x1 = linspace(-x_bar, x_bar, 100);
x2 = linspace(-x_bar, x_bar, 100);
u = linspace(-u_bar,u_bar);
u_max = u_bar;
u_min  = -u_bar;
T = 0.1;
A = [1, T; 0,1];
B = [(T^2) / 2; T];
N = 10;

x_k = combvec(x1, x2);
result = zeros(2,1);
for i = 1:length(x_k)
    x_state = x_k(:,i);
    x_state_change_1 = x_state;
    x_state_change_2 = x_state;
    for j = 2 : N + 1
        x_state_change_1 = A * x_state_change_1 + B * u_max;
        x_state_change_2 = A * x_state_change_2 + B * u_min;
        if ( ~all(x_state_change_1 < [x_bar; x_bar])  || ~all(x_state_change_2 > [-x_bar; -x_bar]) )
                break
        end
        if j == N + 1
            result = [result x_state];
        end
    end
end
len = length(result);
result = result(:, 2:len);
figure
plot(result(1,:), result(2,:), '.')
xlim([-5.0 5.0])
ylim([-5.0 5.0])
%% verifying
x_bar = 5;
u_bar = 1;
u = linspace(-u_bar,u_bar);
u_max = u_bar;
u_min  = -u_bar;
T = 0.1;
A = [1, T; 0,1];
B = [(T^2) / 2; T];
N = 10;

x_k = [-0.05051;3.99];
result = zeros(2,1);
x_state_change_1 = x_k;
x_state_change_2 = x_k;
count = 0;
for j = 2 : N + 1
     x_state_change_1 = A * x_state_change_1 + B * u_max;
     x_state_change_2 = A * x_state_change_2 + B * u_min;
     if ( ~all(x_state_change_1 < [x_bar; x_bar])  || ~all(x_state_change_2 > [-x_bar; -x_bar]) )
             break
     end
     x_state_change_1
     x_state_change_2
     count = count+1
end
%% problem 1.4
close all;
clear all;
x_bar = 5;
u_bar = 1;
u_max = u_bar;
u_min  = -u_bar;

T = 0.1;
A = [1, T; 0,1];
B = [(T^2) / 2; T];

Q = eye(2);
R = 0.1;
S = 10 * eye(2);

N = 10;

xo = [2;0];
x = zeros(2,1,20);
x(:,:,1) = xo;

r = zeros(2,1,20);
r(:,:,1) = xo;
for i = 2:19
    if i <= N
        r(:,:,i) = ((N - (i- 1)) / N) .* xo;
    else
        r(:,:,i) = 0;
    end
end

state_k = xo - r(:,:,1);

for i = 2:N
    
    ref = [0;0;];
    c_bar = [0;0];
    for k = 1:N
        ref = [ref;r(:,:, i+k)];
        c_bar = [c_bar; A^k * r(:,:,i)];
    end
    

    f_bar = [eye(2); A; A^2;A^3;A^4;A^5;A^6;A^7;A^8;A^9;A^10;]*state_k;
    Z = zeros(2,1);
    b_bar = [Z          Z Z Z Z Z Z Z Z Z; 
                  B          Z Z Z Z Z Z Z Z Z;
                  A^1*B          B          Z Z Z Z Z Z Z Z;      
                  A^2*B  A^1*B          B          Z Z Z Z Z Z Z;
                  A^3*B  A^2*B  A^1*B          B          Z Z Z Z Z Z;
                  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z Z Z Z;
                  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z Z Z;
                  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B          Z Z Z;
                  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z;
                  A^8*B  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B  Z;
                  A^9*B  A^8*B  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B  B];

    q_bar = blkdiag(Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,S);
    r_bar = diag([R,R,R,R,R,R,R,R,R,R]);
    Lu1 = eye(N);
    Lu2 = -eye(N);
    Lx1 = eye(2*(N+1));
    Lx2 = -eye(2*(N+1));

    %quadratic programming
    Q_bar = b_bar.' * q_bar* b_bar + r_bar;
    C_T = f_bar.' * q_bar * b_bar - ref.' * q_bar * b_bar + c_bar.' * q_bar * b_bar; 
    C= C_T.';
    A_q = [Lu1;Lu2;Lx1*b_bar;Lx2*b_bar];
    tmp1 = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx1*(f_bar - ref + c_bar);
    tmp2 = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx2*(f_bar - ref + c_bar);
    b_q = [u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               tmp1;
               tmp2
               ];
    [U,~,~,~,~] = quadprog(Q_bar, C,A_q,b_q);
    %system dynamics
    x(:,:,i) = A * x(:,:,i - 1) + B * U(1);
    state_k = x(:,:,i) - r(:,:,i);
    
    %plan trajectories
    plan_trajectory = zeros(2,1,N+1);
    plan_trajectory(:,:,1) = x(:,:,i-1);
    for j = 2:N+1
        plan_trajectory(:,:,j) = A * plan_trajectory(:,:,j - 1) + B * U(j-1);
    end
    figure(1);
    plot( i - 2: i+N - 3, squeeze(plan_trajectory(1,1,1:N)));
    hold on
    figure(2);
    plot( i - 2:i+N - 3, squeeze(plan_trajectory(2,1,1:N)));
    hold on
    

end

figure(1);
plot(0:N - 1, squeeze(x(1,1,1:N)));
figure(2);
plot(0:N - 1, squeeze(x(2,1,1:N)));
figure(1);
plot(0:19, squeeze(r(1,1,1:20)));
legend('MPC','MPC','MPC','MPC','MPC','MPC','MPC','MPC','MPC','executed trajectory','reference trajectory');
figure(2);
plot(0:19, squeeze(r(2,1,1:20)));
legend('MPC','MPC','MPC','MPC','MPC','MPC','MPC','MPC','MPC','executed trajectory','reference trajectory');
%% problem 2.1
clear all
close all
% w = rand(1,10);
w = [0.157613081677548,0.970592781760616,0.957166948242946,0.485375648722841,0.800280468888800,
        0.141886338627215,0.421761282626275,0.915735525189067,0.792207329559554,0.959492426392903];
x_bar = 5;
u_bar = 1;

T = 0.1;
A = [1, T; 0,1];
B = [(T^2) / 2; T];

Q = eye(2);
R = 0.1;
S = 10 * eye(2);

N = 10;

xo = [2;0];
x = zeros(2,1,20);
x(:,:,1) = xo;

r = zeros(2,1,20);
r(:,:,1) = xo;
for i = 2:19
    if i <= N
        r(:,:,i) = ((N - (i- 1)) / N) .* xo;
    else
        r(:,:,i) = 0;
    end
end

state_k = xo - r(:,:,1);
control = zeros(2,1,N-1);
for i = 2:N
    
    ref = [0;0;];
    c_bar = [0;0];
    for k = 1:N
        ref = [ref;r(:,:, i+k)];
        c_bar = [c_bar; A^k * r(:,:,i)];
    end
    

    f_bar = [eye(2); A; A^2;A^3;A^4;A^5;A^6;A^7;A^8;A^9;A^10;]*state_k;
    Z = zeros(2,1);
    b_bar = [Z          Z Z Z Z Z Z Z Z Z; 
                  B          Z Z Z Z Z Z Z Z Z;
                  A^1*B          B          Z Z Z Z Z Z Z Z;      
                  A^2*B  A^1*B          B          Z Z Z Z Z Z Z;
                  A^3*B  A^2*B  A^1*B          B          Z Z Z Z Z Z;
                  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z Z Z Z;
                  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z Z Z;
                  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B          Z Z Z;
                  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z;
                  A^8*B  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B  Z;
                  A^9*B  A^8*B  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B  B];

    q_bar = blkdiag(Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,S);
    r_bar = diag([R,R,R,R,R,R,R,R,R,R]);
    Lu1 = eye(N);
    Lu2 = -eye(N);
    Lx1 = eye(2*(N+1));
    Lx2 = -eye(2*(N+1));

    %quadratic programming
    Q_bar = b_bar.' * q_bar* b_bar + r_bar;
    C_T = f_bar.' * q_bar * b_bar - ref.' * q_bar * b_bar + c_bar.' * q_bar * b_bar; 
    C= C_T.';
    A_q = [Lu1;Lu2;Lx1*b_bar;Lx2*b_bar];
    tmp1 = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx1*(f_bar - ref + c_bar);
    tmp2 = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx2*(f_bar - ref + c_bar);
    b_q = [u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               tmp1;
               tmp2
               ];
    %designed control
    [U,~,~,~,~] = quadprog(Q_bar, C,A_q,b_q);
    control(:,:,i - 1) = U(1);
    
    %system dynamics
    x(:,:,i) = A * x(:,:,i - 1) + B * (U(1) + w(i - 1));
    state_k = x(:,:,i) - r(:,:,i);
    
    %plan trajectories
    plan_trajectory = zeros(2,1,N+1);
    plan_trajectory(:,:,1) = x(:,:,i-1);
    for j = 2:N+1
        plan_trajectory(:,:,j) = A * plan_trajectory(:,:,j - 1) + B * (U(j-1) + w(i - 1));
    end
    figure(1);
    plot( i - 2: i+N - 3, squeeze(plan_trajectory(1,1,1:N)));
    hold on
    figure(2);
    plot( i - 2:i+N - 3, squeeze(plan_trajectory(2,1,1:N)));
    hold on
    

end

figure(1);
plot(0:N - 1, squeeze(x(1,1,1:N)));
legend('MPC','MPC','MPC','MPC','MPC','MPC','MPC','MPC','MPC','executed trajectory');
figure(2);
plot(0:N - 1, squeeze(x(2,1,1:N)));
legend('MPC','MPC','MPC','MPC','MPC','MPC','MPC','MPC','MPC','executed trajectory');
figure(3);
plot(0:8, squeeze(control(1,1,1:9)));
figure(4);
plot(0:8, squeeze(control(2,1,1:9)));

%% problem 2.2
clear all
close all
w = [0.157613081677548,0.970592781760616,0.957166948242946,0.485375648722841,0.800280468888800,
        0.141886338627215,0.421761282626275,0.915735525189067,0.792207329559554,0.959492426392903];
x_bar = 5;
u_bar = 1;

T = 0.1;
A = [1, T; 0,1];
B = [(T^2) / 2; T];

Q = eye(2);
R = 0.1;
S = 10 * eye(2);

N = 10;

xo = [2;0];
x_star = zeros(2,1,10);
x_noisy = zeros(2,1,10);
x_star(:,:,1) = xo;
x_noisy(:,:,1) = xo;

r = zeros(2,1,20);
r(:,:,1) = xo;
for i = 2:19
    if i <= N
        r(:,:,i) = ((N - (i- 1)) / N) .* xo;
    else
        r(:,:,i) = 0;
    end
end


state_k_star = xo - r(:,:,1);
state_k_noisy = xo - r(:,:,1);
control_star = zeros(2,1,N-1);
for i = 2:N
    
    ref = [0;0;];
    c_bar = [0;0];
    for k = 1:N
        ref = [ref;r(:,:, i+k)];
        c_bar = [c_bar; A^k * r(:,:,i)];
    end
    

    f_bar_star = [eye(2); A; A^2;A^3;A^4;A^5;A^6;A^7;A^8;A^9;A^10;]*state_k_star;
    f_bar_noisy = [eye(2); A; A^2;A^3;A^4;A^5;A^6;A^7;A^8;A^9;A^10;]*state_k_noisy;
    Z = zeros(2,1);
    b_bar = [Z          Z Z Z Z Z Z Z Z Z; 
                  B          Z Z Z Z Z Z Z Z Z;
                  A^1*B          B          Z Z Z Z Z Z Z Z;      
                  A^2*B  A^1*B          B          Z Z Z Z Z Z Z;
                  A^3*B  A^2*B  A^1*B          B          Z Z Z Z Z Z;
                  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z Z Z Z;
                  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z Z Z;
                  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B          Z Z Z;
                  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z;
                  A^8*B  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B  Z;
                  A^9*B  A^8*B  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B  B];

    q_bar = blkdiag(Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,S);
    r_bar = diag([R,R,R,R,R,R,R,R,R,R]);
    Lu1 = eye(N);
    Lu2 = -eye(N);
    Lx1 = eye(2*(N+1));
    Lx2 = -eye(2*(N+1));

    %quadratic programming
    Q_bar = b_bar.' * q_bar* b_bar + r_bar;
    C_T_star = f_bar_star.' * q_bar * b_bar - ref.' * q_bar * b_bar + c_bar.' * q_bar * b_bar; 
    C_star= C_T_star.';
    C_T_noisy = f_bar_noisy.' * q_bar * b_bar - ref.' * q_bar * b_bar + c_bar.' * q_bar * b_bar; 
    C_noisy= C_T_noisy.';
    
    A_q = [Lu1;Lu2;Lx1*b_bar;Lx2*b_bar];

    tmp1_star = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx1*(f_bar_star - ref + c_bar);
    tmp2_star = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx2*(f_bar_star - ref + c_bar);
    b_q_star = [u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               tmp1_star;
               tmp2_star
               ];
    tmp1_noisy = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx1*(f_bar_noisy - ref + c_bar);
    tmp2_noisy = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx2*(f_bar_noisy - ref + c_bar);
    b_q_noisy = [u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               tmp1_noisy;
               tmp2_noisy
               ];
           
    %designed control
    [U_star,~,~,~,~] = quadprog(Q_bar, C_star,A_q,b_q_star);
    control_star(:,:,i - 1) = U_star(1);
    [U_noisy,~,~,~,~] = quadprog(Q_bar, C_noisy,A_q,b_q_noisy);
    
    %system dynamics
    x_star(:,:,i) = A * x_star(:,:,i - 1) + B * U_star(1);
    state_k_star = x_star(:,:,i) - r(:,:,i);
    
    x_noisy(:,:,i) = A * x_noisy(:,:,i - 1) + B * (U_noisy(1) + w(i - 1));
    state_k_noisy = x_noisy(:,:,i) - r(:,:,i);

end

%plot states trajectories
figure(1);
title('first state')
plot(0:N - 1, squeeze(x_star(1,1,1:N)));
hold on;
plot(0:N - 1, squeeze(x_noisy(1,1,1:N)));
figure(2);
title('second state')
plot(0:N - 1, squeeze(x_star(2,1,1:N)));
hold on;
plot(0:N - 1,squeeze(x_noisy(2,1,1:N)));

%plot errors
error = zeros(2,1,N);
error(:,:,:) = x_noisy(:,:,:) - x_star(:,:,:);
figure(3);
title('errors')
plot(0:N - 1, squeeze(error(1,1,1:N)));
hold on
plot(0:N - 1, squeeze(error(2,1,1:N)));
legend("state_{11}","state_{21}")

%% problem 2.3
% clear all
% close all
% w = [0.157613081677548,0.970592781760616,0.957166948242946,0.485375648722841,0.800280468888800,
%         0.141886338627215,0.421761282626275,0.915735525189067,0.792207329559554,0.959492426392903];
% x_bar = 5;
% u_bar = 1;
% 
% T = 0.1;
% A = [1, T; 0,1];
% B = [(T^2) / 2; T];
% 
% Q = eye(2);
% R = 0.1;
% S = 10 * eye(2);
% 
% N = 10;
% 
% xo = [2;0];
% x_star = zeros(2,1,10);
% x_noisy = zeros(2,1,10);
% x_star(:,:,1) = xo;
% x_noisy(:,:,1) = xo;
% 
% r = zeros(2,1,20);
% r(:,:,1) = xo;
% for i = 2:19
%     if i <= N
%         r(:,:,i) = ((N - (i- 1)) / N) .* xo;
%     else
%         r(:,:,i) = 0;
%     end
% end
% 
% 
% state_k_star = xo - r(:,:,1);
% state_k_noisy = xo - r(:,:,1);
% control_star = zeros(1,1,N-1);
% control_noisy = zeros(1,1,N - 1);
% for i = 2:N
%     
%     ref = [0;0;];
%     c_bar = [0;0];
%     for k = 1:N
%         ref = [ref;r(:,:, i+k)];
%         c_bar = [c_bar; A^k * r(:,:,i)];
%     end
%     
% 
%     f_bar_star = [eye(2); A; A^2;A^3;A^4;A^5;A^6;A^7;A^8;A^9;A^10;]*state_k_star;
%     f_bar_noisy = [eye(2); A; A^2;A^3;A^4;A^5;A^6;A^7;A^8;A^9;A^10;]*state_k_noisy;
%     Z = zeros(2,1);
%     b_bar = [Z          Z Z Z Z Z Z Z Z Z; 
%                   B          Z Z Z Z Z Z Z Z Z;
%                   A^1*B          B          Z Z Z Z Z Z Z Z;      
%                   A^2*B  A^1*B          B          Z Z Z Z Z Z Z;
%                   A^3*B  A^2*B  A^1*B          B          Z Z Z Z Z Z;
%                   A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z Z Z Z;
%                   A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z Z Z;
%                   A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B          Z Z Z;
%                   A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z;
%                   A^8*B  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B  Z;
%                   A^9*B  A^8*B  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B  B];
% 
%     q_bar = blkdiag(Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,S);
%     r_bar = diag([R,R,R,R,R,R,R,R,R,R]);
%     Lu1 = eye(N);
%     Lu2 = -eye(N);
%     Lx1 = eye(2*(N+1));
%     Lx2 = -eye(2*(N+1));
% 
%     %quadratic programming
%     Q_bar = b_bar.' * q_bar* b_bar + r_bar;
%     C_T_star = f_bar_star.' * q_bar * b_bar - ref.' * q_bar * b_bar + c_bar.' * q_bar * b_bar; 
%     C_star= C_T_star.';
%     C_T_noisy = f_bar_noisy.' * q_bar * b_bar - ref.' * q_bar * b_bar + c_bar.' * q_bar * b_bar; 
%     C_noisy= C_T_noisy.';
%     
%     A_q = [Lu1;Lu2;Lx1*b_bar;Lx2*b_bar];
%     tmp1_star = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx1*f_bar_star;
%     tmp2_star = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx2*f_bar_star;
%     b_q_star = [u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
%                u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
%                tmp1_star;
%                tmp2_star
%                ];
%     tmp1_noisy = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx1*f_bar_noisy;
%     tmp2_noisy = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx2*f_bar_noisy;
%     b_q_noisy = [u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
%                u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
%                tmp1_noisy;
%                tmp2_noisy
%                ];
%            
%     %designed control
%     [U_star,~,~,~,~] = quadprog(Q_bar, C_star,A_q,b_q_star);
%     control_star(:,:,i - 1) = U_star(1);
%     [U_noisy,~,~,~,~] = quadprog(Q_bar, C_noisy,A_q,b_q_noisy);
%     control_noisy(:,:,i - 1) = U_noisy(1);
%     
%     %system dynamics
%     x_star(:,:,i) = A * x_star(:,:,i - 1) + B * U_star(1);
%     state_k_star = x_star(:,:,i) - r(:,:,i);
%     
%     
%     x_noisy(:,:,i) = A * x_noisy(:,:,i - 1) + B * (U_noisy(1) + w(i - 1));
%     state_k_noisy = x_noisy(:,:,i) - r(:,:,i);
%     
% 
% end
% 
% %plot errors
% error = zeros(2,1,N);
% error(:,:,:) = x_noisy(:,:,:) - x_star(:,:,:);
% figure(3);
% title('errors')
% plot(0:N - 1, squeeze(error(1,1,1:N)));
% hold on
% plot(0:N - 1, squeeze(error(2,1,1:N)));
% hold on
% 
% 
% %second iteration to fix the error
% state_k_noisy = xo - r(:,:,1);
% L = zeros(2,1,N);
% for i = 2:N
%     %find optimal L
%     tmp1 = x_star(:,:,i);
%     tmp2 = x_noisy(:,:,i - 1);
%     err = error(:,:,i - 1);
%     u = control_noisy(:,:,i - 1);
%     L_guess = [-1  -1];
%     fun = @(L_itr) cal(L_itr, tmp1, tmp2, A, B, u, err,w(i - 1));
%     [L_min, fval, exitflag]= fminsearch(fun, [-1 -1]);
%     L(:,:,i - 1) = L_min(1);
%     L(:,:,i - 1) = L_min(2);
%     
%     %system dynamics
%     new_control_law = control_noisy(i - 1) + L (:,:,i - 1).' * error(:,:,i - 1);
%     x_noisy(:,:,i) = A * x_noisy(:,:,i - 1) + B * (new_control_law + w(i - 1));
%     state_k_noisy = x_noisy(:,:,i) - r(:,:,i);
% 
% end
% 
% %plot errors
% error = zeros(2,1,N);
% error(:,:,:) = x_noisy(:,:,:) - x_star(:,:,:);
% figure(3);
% title('errors')
% plot(0:N - 1, squeeze(error(1,1,1:N)));
% hold on
% plot(0:N - 1, squeeze(error(2,1,1:N)));
% legend("state_{11}","state_{21}","new state_{11}","new state_{21}")
% 
% function f = cal(Li, tmp1, tmp2, A, B, u, err, noise)
%      tmp = A * tmp2 + B * (u + Li * err  + noise) - tmp1;
%      f = abs(tmp(1)) + abs(tmp(2));
% end


%% problem 2.4
clear all
close all
w = [0.157613081677548,0.970592781760616,0.957166948242946,0.485375648722841,0.800280468888800,
        0.141886338627215,0.421761282626275,0.915735525189067,0.792207329559554,0.959492426392903];
x_bar = 5;
u_bar = 1;

T = 0.1;
A = [1, T; 0,1];
B = [(T^2) / 2; T];

Q = eye(2);
R = 0.1;
S = 10 * eye(2);

N = 10;

xo = [2;0];
x_star = zeros(2,1,10);
x_noisy = zeros(2,1,10);
x_star(:,:,1) = xo;
x_noisy(:,:,1) = xo;

r = zeros(2,1,20);
r(:,:,1) = xo;
for i = 2:19
    if i <= N
        r(:,:,i) = ((N - (i- 1)) / N) .* xo;
    else
        r(:,:,i) = 0;
    end
end


state_k_star = xo - r(:,:,1);
state_k_noisy = xo - r(:,:,1);
control_star = zeros(1,1,N-1);
control_noisy = zeros(1,1,N - 1);
for i = 2:N
    
    ref = [0;0;];
    c_bar = [0;0];
    for k = 1:N
        ref = [ref;r(:,:, i+k)];
        c_bar = [c_bar; A^k * r(:,:,i)];
    end
    

    f_bar_star = [eye(2); A; A^2;A^3;A^4;A^5;A^6;A^7;A^8;A^9;A^10;]*state_k_star;
    f_bar_noisy = [eye(2); A; A^2;A^3;A^4;A^5;A^6;A^7;A^8;A^9;A^10;]*state_k_noisy;
    Z = zeros(2,1);
    b_bar = [Z          Z Z Z Z Z Z Z Z Z; 
                  B          Z Z Z Z Z Z Z Z Z;
                  A^1*B          B          Z Z Z Z Z Z Z Z;      
                  A^2*B  A^1*B          B          Z Z Z Z Z Z Z;
                  A^3*B  A^2*B  A^1*B          B          Z Z Z Z Z Z;
                  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z Z Z Z;
                  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z Z Z;
                  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B          Z Z Z;
                  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B         Z Z;
                  A^8*B  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B          B  Z;
                  A^9*B  A^8*B  A^7*B  A^6*B  A^5*B  A^4*B  A^3*B  A^2*B  A^1*B  B];

    q_bar = blkdiag(Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,S);
    r_bar = diag([R,R,R,R,R,R,R,R,R,R]);
    Lu1 = eye(N);
    Lu2 = -eye(N);
    Lx1 = eye(2*(N+1));
    Lx2 = -eye(2*(N+1));

    %quadratic programming
    Q_bar = b_bar.' * q_bar* b_bar + r_bar;
    C_T_star = f_bar_star.' * q_bar * b_bar - ref.' * q_bar * b_bar + c_bar.' * q_bar * b_bar; 
    C_star= C_T_star.';
    C_T_noisy = f_bar_noisy.' * q_bar * b_bar - ref.' * q_bar * b_bar + c_bar.' * q_bar * b_bar; 
    C_noisy= C_T_noisy.';
    
    A_q = [Lu1;Lu2;Lx1*b_bar;Lx2*b_bar];
    tmp1_star = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx1*f_bar_star;
    tmp2_star = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx2*f_bar_star;
    b_q_star = [u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               tmp1_star;
               tmp2_star
               ];
    tmp1_noisy = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx1*f_bar_noisy;
    tmp2_noisy = [x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar;x_bar] - Lx2*f_bar_noisy;
    b_q_noisy = [u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;u_bar;
               tmp1_noisy;
               tmp2_noisy
               ];
           
    %designed control
    [U_star,~,~,~,~] = quadprog(Q_bar, C_star,A_q,b_q_star);
    control_star(:,:,i - 1) = U_star(1);
    [U_noisy,~,~,~,~] = quadprog(Q_bar, C_noisy,A_q,b_q_noisy);
    control_noisy(:,:,i - 1) = U_noisy(1);
    
    %system dynamics
    x_star(:,:,i) = A * x_star(:,:,i - 1) + B * U_star(1);
    state_k_star = x_star(:,:,i) - r(:,:,i);
    
    
    x_noisy(:,:,i) = A * x_noisy(:,:,i - 1) + B * (U_noisy(1) + w(i - 1));
    state_k_noisy = x_noisy(:,:,i) - r(:,:,i);
    

end

%plot errors
error = zeros(2,1,N);
error(:,:,:) = x_noisy(:,:,:) - x_star(:,:,:);
figure(3);
title('errors')
plot(0:N - 1, squeeze(error(1,1,1:N)));
hold on
plot(0:N - 1, squeeze(error(2,1,1:N)));
hold on


%second iteration to fix the error
state_k_noisy = xo - r(:,:,1);
L = zeros(2,1,N);
iteration = 10;
Q = 0.99;   % bad at 0.8 + 
for itr = 1 : iteration
    for i = 2:N
        %find optimal L
        tmp1 = x_star(:,:,i);
        tmp2 = x_noisy(:,:,i - 1);
        err = error(:,:,i - 1);
        u = control_noisy(:,:,i - 1);
        L_guess = [-1  -1];
        fun = @(L_itr) cal(L_itr, tmp1, tmp2, A, B, u, err,w(i - 1));
        [L_min, fval, exitflag]= fminsearch(fun, [-1 -1]);
        L(:,:,i - 1) = L_min(1);
        L(:,:,i - 1) = L_min(2);

        %system dynamics
        new_control_law = control_noisy(i - 1) + L (:,:,i - 1).' * error(:,:,i - 1);
        x_noisy(:,:,i) = A * x_noisy(:,:,i - 1) + B * (Q* new_control_law + w(i - 1));
        control_noisy(:,:,i - 1) = Q * new_control_law;
        state_k_noisy = x_noisy(:,:,i) - r(:,:,i);
    end
    %plot errors
    error = zeros(2,1,N);
    error(:,:,:) = x_noisy(:,:,:) - x_star(:,:,:);
    figure(3);
    title('errors_{11}')
    plot(0:N - 1, squeeze(error(1,1,1:N)));
    hold on

    figure(4);
    title('errors_{21}')
    plot(0:N - 1, squeeze(error(2,1,1:N)));
    hold on
end





function f = cal(Li, tmp1, tmp2, A, B, u, err, noise)
     tmp = A * tmp2 + B * (u + Li * err  + noise) - tmp1;
     f = abs(tmp(1)) + abs(tmp(2));
end

