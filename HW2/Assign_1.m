% %% problem 2
% velocity = 0
% theta = 0* pi / 180;
% 
%  A =[0, 0, cos(theta), -velocity * sin(theta); 0, 0, sin(theta), velocity * cos(theta); 0, 0, 0, 0; 0, 0, 0, 0];
%  B = [0, 0; 0, 0; 1, 0; 0, 1];
%  Q = eye(4);
%  R = 0.1 * eye(2);
%  T = 10;
% 
% x0 = [0;0;0;0];
% [K, P] = lqr(A, B, Q, R);
% global control_gain
% control_gain = K
% %% ODE problem 2
% [t, x] = ode45(@rhs, [0,10],x0);
% figure;
% plot(x(:,1),'k')
% hold on;
% plot(x(:,2), 'b')
% hold on;
% plot(x(:,3), 'r')
% hold on;
% plot(x(:,4), 'g')
% legend("x_1", "x_2","x_3","x_4")
% function xdot = rhs(t,x)
%     global control_gain
%     tmp = - control_gain * (x - [1;1;0;0])
%     xdot_1 = x(3) * cos(x(4));
%     xdot_2 = x(3) * sin(x(4));
%     xdot_3 = tmp(1);
%     xdot_4 = tmp(2);
%     xdot = [xdot_1; xdot_2; xdot_3; xdot_4];
% end
   
%% problem 3.5
syms theta velocity x1 x2 x3 x4 xn1 xn2 xn3 xn4 n1


f_x = [1, 0, 0.1 * cos(x4), -0.1* x3 * sin(x4); 0 , 1, 0.1 * sin(x4), 0.1* x3 * cos(x4); 0,0,1,0; 0,0,0,1];
Q = eye(4);
R = 0.1 * eye(2);
S = 10 * eye(4);
N = 100;
xg = [10;10; 0; pi/2];
% xg = [xg1; xg2; xg3; xg4];



state_N_1 = [x1;x2;x3;x4];
state_N = [xn1; xn2; xn3; xn4];
% state_N = [x1 + 0.1 * x3 * cos(x4); x2 + 0.1 * x3 * sin(x4); x3 + 0.1 * u1; x4 + 0.1 * u2];

co_state = sym(zeros(4,1,N));
co_state(:,:,N) = S * (state_N -  xg);
c1 = sym(zeros(1,N)); c1(N) = co_state(1,1,N);
c2 = sym(zeros(1,N)); c2(N) = co_state(2,1,N);
c3 = sym(zeros(1,N)); c3(N) = co_state(3,1,N);
c4 = sym(zeros(1,N)); c4(N) = co_state(4,1,N);

u = -inv(R) * [0,0,0.1,0; 0,0,0,0.1] * co_state(:,:,N);
% plug in co_state(:,:,N) and represent u in terms of co_state(:,:, N - 1)

co_state(:,:,N - 1) = Q * (state_N_1 - xg) +  f_x.' * co_state(:,:,N);
co_state(:,:,N - 1) = subs(co_state(:,:,N - 1), [xn1, xn2, xn3, xn4], [x1 + 0.1 * x3 * cos(x4), x2 + 0.1 * x3 * sin(x4), x3 + 0.1 * (-5 * x3), x4 + 0.1 * (-5 * x4)]);
co_state(:,:,N - 1) = taylor(co_state(:,1,N - 1), [x1;x2;x3;x4],'ExpansionPoint', xg,'Order',2);

%%

B = [0, 0; 0, 0; 0.1, 0; 0, 0.1];
Q = eye(4);
R = 0.1 * eye(2);
S = 10 * eye(4);
N = 100;

x0 = [0;0;0;0];
xg = [10;10; 0; pi/2];

state = zeros(4,1,N);
state(:,:,1) = x0;
s11 = zeros(1,N); s11(1) = state(1,1,1);
s21 = zeros(1,N); s21(1) = state(2,1,1);
s31 = zeros(1,N); s31(1) = state(3,1,1);
s41 = zeros(1,N); s41(1) = state(4,1,1);


co_state = zeros(4,1,N);   
c1 = zeros(1,N); c1(1) = 11 * s11(1) - 110;
c2 = zeros(1,N);c2(1) = 11 * s21(1) + s31(1) - 110;
c3 = zeros(1,N);c3(1) = s21(1) + (61 * s31(1) / 10) - 10;
c4 = zeros(1,N);c4(1) = 6* s41(1) - (11 * pi) / 2;

u = zeros(2,1,N);
u(:,:,1) = -inv(R) * [0,0,0,0.1; 0,0,0.1,0] * co_state(:,:,2);
u1 = zeros(1,N); u1(1) = u(1,1,1);
u2 = zeros(1,N); u2(1) = u(2,1,1);
 for i = 2:N
    state(:,:,i) = [0.1 * s31(i - 1) * cos(s41(i - 1)); 0.1* s31(i - 1) * sin(s41(i - 1)); 0 ; 0] + state(:,:,i - 1)+ B * u(:,:,i - 1);
    s11(i) = state(1,1,i);
    s21(i) = state(1,1,i);
    s31(i) = state(3,1,i);
    s41(i) = state(4,1,i);
        
    c1(i) = 11 * s11(i) - 110;
    c2(i) = 11 * s21(i) + s31(i) - 110;
    c3(i) = s21(i) + (61 * s31(i) / 10) - 10;
    c4(i) = 6* s41(i) - (11 * pi) / 2;
    co_state(:,:,i) = [c1(i);c2(i);c3(i);c4(i)];
    
    u(:,:,i) = -inv(R) * [0,0,0.1,0; 0,0,0,0.1] * co_state(:,:,i);
    u1(i) = u(1,1,i);
    u2(i) = u(2,1,i);
 end


figure;
hold on;
plot(s11)
plot(s21)
plot(s31)
plot(s41)
legend("s_{1}", "s_{2}", "s_{3}", "s_{4}")

figure;
hold on;
plot(c1)
plot(c2)
plot(c3)
plot(c4)
legend("c_{1}", "c_{2}", "c_{3}", "c_{4}")
figure;
hold on;
plot(u1)
plot(u2)
legend("u_{1}", "u_{2}")


%% problem 3.6
% clc 
% clear all
% close all
% velocity = 0;
% theta = 90 * pi / 180;
%  
% A =[1, 0, 0.1* cos(theta), 0.1 * -velocity * sin(theta); 0, 1, 0.1 * sin(theta), 0.1 * velocity * cos(theta); 0, 0, 1, 0; 0, 0, 0, 1];
% B = [0, 0; 0, 0; 0.1, 0; 0, 0.1];
% Q = eye(4);
% R = 0.1 * eye(2);
% S = 10 * eye(4);
% N = 1000;
% 
% P = cell(N);
% P{N} = S;
% p11 = zeros(1,N); p11(N) = P{N}(1,1);
% p12 = zeros(1,N); p12(N) = P{N}(1,2);
% p13 = zeros(1,N); p13(N) = P{N}(1,3);
% p14 = zeros(1,N); p14(N) = P{N}(1,4);
% p21 = zeros(1,N); p21(N) = P{N}(2,1);
% p22 = zeros(1,N); p22(N) = P{N}(2,2);
% p23 = zeros(1,N); p23(N) = P{N}(2,3);
% p24 = zeros(1,N); p24(N) = P{N}(2,4);
% p31 = zeros(1,N); p31(N) = P{N}(3,1);
% p32 = zeros(1,N); p32(N) = P{N}(3,2);
% p33 = zeros(1,N); p33(N) = P{N}(3,3);
% p34 = zeros(1,N); p34(N) = P{N}(3,4);
% p41 = zeros(1,N); p41(N) = P{N}(4,1);
% p42 = zeros(1,N); p42(N) = P{N}(4,2);
% p43 = zeros(1,N); p43(N) = P{N}(4,3);
% p44 = zeros(1,N); p44(N) = P{N}(4,4);
% 
% % Solving Riccati Equation
% for i = N-1:-1:1
%     P{i} = Q + A' * P{i+1} * A - A' * P{i+1} * B * inv(R + B' * P{i+1} * B) * B' * P{i+1} * A;
%     p11(i) = P{i}(1,1);
%     p12(i) = P{i}(1,2);
%     p13(i) = P{i}(1,3);
%     p14(i) = P{i}(1,4);
%     p21(i) = P{i}(2,1);
%     p22(i) = P{i}(2,2);
%     p23(i) = P{i}(2,3);
%     p24(i) = P{i}(2,4);    
%     p31(i) = P{i}(3,1);
%     p32(i) = P{i}(3,2);
%     p33(i) = P{i}(3,3);    
%     p34(i) = P{i}(3,4);
%     p41(i) = P{i}(4,1);
%     p42(i) = P{i}(4,2);
%     p43(i) = P{i}(4,3);
%     p44(i) = P{i}(4,4);
%     
% end
% 
% % % Plot
% figure; hold on;
% plot(p11)
% plot(p12)
% plot(p21)
% plot(p22)
% legend("p_{11}", "p_{12}", "p_{21}", "p_{22}")
% 
% x0 = [0;0;0;0];
% xg = [10; 10; 0 ; pi / 2];
% 
% state = cell(N);
% state{1} = x0;
% s11 = zeros(1,N); s11(1) = state{1}(1,1);
% s21 = zeros(1,N); s21(1) = state{1}(2,1);
% s31 = zeros(1,N); s31(1) = state{1}(3,1);
% s41 = zeros(1,N); s41(1) = state{1}(4,1);
% 
% 
% co_state = cell(N);
% co_state{1} = P{1} * state{1};
% c1 = zeros(1,N); c1(1) = co_state{1}(1,1);
% c2 = zeros(1,N); c2(1) = co_state{1}(2,1);
% c3 = zeros(1,N); c3(1) = co_state{1}(3,1);
% c4 = zeros(1,N); c4(1) = co_state{1}(4,1);
% 
% K = cell(N);
% K{1} = - inv(B' * P{1} * B + R) * B' * P{1} * A;
% u = cell(N);
% u{1} = K{1} * (state{1} - xg);
% u1 = zeros(1,N); u1(1) = u{1}(1,1);
% u2 = zeros(1,N); u2(1) = u{1}(2,1);
% 
% 
% for i = 2:N-1
%     state{i} = A * state{i - 1} + B * u{i - 1};
%     s11(i) = state{i}(1,1);
%     s21(i) = state{i}(2,1);
%     s31(i) = state{i}(3,1);
%     s41(i) = state{i}(4,1);
%     
%     co_state{i} = P{i} * state{i};
%     c1(i) = co_state{i}(1,1);
%     c2(i) = co_state{i}(2,1);
%     c3(i) = co_state{i}(3,1);
%     c4(i) = co_state{i}(4,1);
%     
%     
%     K{i} =  - inv(B' * P{i+1} * B + R) * B' * P{i+1} * A;
% 
%     u{i} = K{i} * (state{i} - xg);
%     u1(i) = u{i}(1,1);
%     u2(i) = u{i}(2,1);
% 
% end
% state{N} = A * state{N - 1} + B * u{N - 1};
% s11(N) = state{N}(1,1);
% s21(N) = state{N}(2,1);
% s31(N) = state{N}(3,1);
% s41(N) = state{N}(4,1);
% co_state{N} = P{N} * state{N};
% c1(N) = co_state{N}(1,1);
% c2(N) = co_state{N}(2,1);
% c3(N) = co_state{N}(3,1);
% c4(N) = co_state{N}(4,1);
% figure;
% hold on;
% plot(s11)
% plot(s21)
% plot(s31)
% plot(s41)
% legend("s_{1}", "s_{2}", "s_{3}", "s_{4}")
% 
% figure;
% hold on;
% plot(c1)
% plot(c2)
% plot(c3)
% plot(c4)
% legend("c_{1}", "c_{2}", "c_{3}", "c_{4}")
% figure;
% hold on;
% plot(u1)
% plot(u2)
% legend("u_{1}", "u_{2}")

%% problem 3.7
% velocity = 0;
% theta = 90 * pi / 180;
%  
% A =[1, 0, 0.1* cos(theta), 0.1 * -velocity * sin(theta); 0, 1, 0.1 * sin(theta), 0.1 * velocity * cos(theta); 0, 0, 1, 0; 0, 0, 0, 1];
% B = [0, 0; 0, 0; 0.1, 0; 0, 0.1];
% Q = eye(4);
% R = 0.1 * eye(2);
% S = 10 * eye(4);
% N = 100;
% 
% P = cell(N);
% P{N} = S;
% p11 = zeros(1,N); p11(N) = P{N}(1,1);
% p12 = zeros(1,N); p12(N) = P{N}(1,2);
% p13 = zeros(1,N); p13(N) = P{N}(1,3);
% p14 = zeros(1,N); p14(N) = P{N}(1,4);
% p21 = zeros(1,N); p21(N) = P{N}(2,1);
% p22 = zeros(1,N); p22(N) = P{N}(2,2);
% p23 = zeros(1,N); p23(N) = P{N}(2,3);
% p24 = zeros(1,N); p24(N) = P{N}(2,4);
% p31 = zeros(1,N); p31(N) = P{N}(3,1);
% p32 = zeros(1,N); p32(N) = P{N}(3,2);
% p33 = zeros(1,N); p33(N) = P{N}(3,3);
% p34 = zeros(1,N); p34(N) = P{N}(3,4);
% p41 = zeros(1,N); p41(N) = P{N}(4,1);
% p42 = zeros(1,N); p42(N) = P{N}(4,2);
% p43 = zeros(1,N); p43(N) = P{N}(4,3);
% p44 = zeros(1,N); p44(N) = P{N}(4,4);
% 
% 
% % Solving Riccati Equation
% for i = N-1:-1:1
%     P{i} = Q + A' * P{i+1} * A - A' * P{i+1} * B * inv(R + B' * P{i+1} * B) * B' * P{i+1} * A;
%     p11(i) = P{i}(1,1);
%     p12(i) = P{i}(1,2);
%     p13(i) = P{i}(1,3);
%     p14(i) = P{i}(1,4);
%     p21(i) = P{i}(2,1);
%     p22(i) = P{i}(2,2);
%     p23(i) = P{i}(2,3);
%     p24(i) = P{i}(2,4);    
%     p31(i) = P{i}(3,1);
%     p32(i) = P{i}(3,2);
%     p33(i) = P{i}(3,3);    
%     p34(i) = P{i}(3,4);
%     p41(i) = P{i}(4,1);
%     p42(i) = P{i}(4,2);
%     p43(i) = P{i}(4,3);
%     p44(i) = P{i}(4,4);
%     
% end
% 
% 
% x0 = [0;0;0;0];
% xg = [10; 10; 0 ; pi / 2];
% 
% state = cell(N);
% state{1} = x0;
% s11 = zeros(1,N); s11(1) = state{1}(1,1);
% s21 = zeros(1,N); s21(1) = state{1}(2,1);
% s31 = zeros(1,N); s31(1) = state{1}(3,1);
% s41 = zeros(1,N); s41(1) = state{1}(4,1);
% 
% 
% co_state = cell(N);
% co_state{1} = P{1} * state{1};
% c1 = zeros(1,N); c1(1) = co_state{1}(1,1);
% c2 = zeros(1,N); c2(1) = co_state{1}(2,1);
% c3 = zeros(1,N); c3(1) = co_state{1}(3,1);
% c4 = zeros(1,N); c4(1) = co_state{1}(4,1);
% 
% K = cell(N);
% K{1} = - inv(B' * P{1} * B + R) * B' * P{1} * A;
% 
% u = cell(N);
% u{1} = K{1} * (state{1} - xg);
% u1 = zeros(1,N); u1(1) = u{1}(1,1);
% u2 = zeros(1,N); u2(1) = u{1}(2,1);
% 
% 
% 
% for i = 2:N - 1
%     state{i} = [0.1 * s31(i - 1) * cos(s41(i - 1)); 0.1* s31(i - 1) * sin(s41(i - 1)); 0 ; 0] + state{i - 1}+ B * u{i - 1};
%     s11(i) = state{i}(1,1);
%     s21(i) = state{i}(2,1);
%     s31(i) = state{i}(3,1);
%     s41(i) = state{i}(4,1);
%     
%     K{i} =  - inv(B' * P{i + 1} * B + R) * B' * P{i + 1} * A;
%     
%     u{i} = K{i} * (state{i} - xg);
%     u1(i) = u{i}(1,1);
%     u2(i) = u{i}(2,1);
% 
% end
% 
% state{N} = [0.1 * s31(N - 1) * cos(s41(N - 1)); 0.1* s31(N - 1) * sin(s41(N - 1)); 0 ; 0] + state{N - 1}+ B * u{N - 1};
% s11(N) = state{N}(1,1);
% s21(N) = state{N}(2,1);
% s31(N) = state{N}(3,1);
% s41(N) = state{N}(4,1);
% 
% figure;
% hold on;
% plot(s11)
% plot(s21)
% plot(s31)
% plot(s41)
% legend("s_{1}", "s_{2}", "s_{3}", "s_{4}")
% 
% 
% figure;
% hold on;
% plot(u1)
% plot(u2)
% legend("u_{1}", "u_{2}")

%% problem 4
% % 
% % 
% A =[1 1; 0 1];
% B = [0.5; 1];
% Q = eye(2);
% N = [1; 1];
% R = 1;
% 
% T = 100;
% 
% P = cell(1);
% P{T} = [1 3;2 1];
% p11 = zeros(1,T); p11(T) = P{T}(1,1);
% p12 = zeros(1,T); p12(T) = P{T}(1,2);
% p21 = zeros(1,T); p21(T) = P{T}(2,1);
% p22 = zeros(1,T); p22(T) = P{T}(2,2);
% 
% for i = T-1:-1:1
%     P{i} =  A.'* P{i+1} *A + 2 * Q + (A.' * P{i+1} * B + 2 * N) * (-inv(2 * R + B.' * P{i+1} * B) * (2 * N.' + B.' * P{i+1} * A)) ;
%     p11(i) = P{i}(1,1);
%     p12(i) = P{i}(1,2);
%     p21(i) = P{i}(2,1);
%     p22(i) = P{i}(2,2);
% end
% control_gain = -inv(2 * R + B.' * P{1} * B) *(2 * N.' + B.' * P{1} * A);
% 
% figure; hold on;
% plot(p11)
% plot(p12)
% plot(p21)
% plot(p22)
% legend("p_{11}", "p_{12}", "p_{21}", "p_{22}")
% [k, p] = dlqr(A, B, 2*Q, 2*R,N);