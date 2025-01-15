clc
clear
close all
%%
syms q1 q2
%Constants
m = 10;
l = 1;
g = 9.81;

m1 = m/2;
m2 = m/2;
%%
%Matricies
M = zeros(2,2);
K = zeros(2,2);
Q = zeros(2,1);

M = l*[m1+m2, m2;...
       m2   , m2];
K = g*[(m1+m2), 0;...
       0      , m2];

%%
% Given mass and stiffness matrices
% M = [10, 5; 5, 5];
% K = [98.1, 0; 0, 49];


% Solve the eigenvalue problem
[V, D] = eig(inv(M) * K);

% Extract the diagonal matrix D containing eigenvalues
eigenvalues = sqrt(diag(D));

% Sort eigenvalues and corresponding eigenvectors in ascending order
[sorted_eigenvalues, index] = sort(eigenvalues);
sorted_eigenvectors = V(:, index);

% Display the sorted eigenvalues (natural frequencies)
disp('Natural Frequencies:');
disp(sorted_eigenvalues);

% Display the sorted eigenvector matrix (mode shapes)
disp('Mode Shapes:');
disp(sorted_eigenvectors);

% Display mode shapes and corresponding natural frequencies
disp('Mode Shapes and Corresponding Natural Frequencies:');
mode1 = sorted_eigenvectors(:, 1);
mode2 = sorted_eigenvectors(:, 2);
w1 = sorted_eigenvalues(1);
w2 = sorted_eigenvalues(2);
for i = 1:length(sorted_eigenvalues)
    natural_frequency = sorted_eigenvalues(i);
    mode_shape = sorted_eigenvectors(:, i);

    fprintf('Mode %d: Natural Frequency = %.4f rad/s, Mode Shape = [%f, %f]\n', i, natural_frequency, mode_shape(1), mode_shape(2));
end

%%
q0 = [ -0.86 ; 0.55];
dq0 = [ -0.81 ; -0.15];


%Find a
n1 = mode1'*M*q0;
d1 = mode1'*M*mode1;
% t1 = mode1'*K*mode1;
% test1 = sqrt(t1/d1)

n2 = mode2'*M*q0;
d2 = mode2'*M*mode2;
% t2 = mode2'*K*mode2;
% test2 = sqrt(t2/d2)

a1 = n1/d1;
a2 = n2/d2;

%Find b

n1 = mode1'*M*dq0;
d1 = mode1'*M*mode1;

n2 = mode2'*M*dq0;
d2 = mode2'*M*mode2;

b1 = n1/(d1*w1);
b2 = n2/(d2*w2);

%Response:
t = 0:0.01:10;
q = a1*mode1*cos(w1*t) + a2*mode2*cos(w2*t) + b1*mode1*sin(w1*t) + b2*mode2*sin(w2*t);
q1_response = q(1,:);
q2_response = q(2,:);

%Double Plots
figure;

subplot(2, 1, 1);
plot(t, rad2deg(q1_response), 'LineWidth', 2, 'DisplayName', '\theta_1');
xlabel('Time [s]');
ylabel('\theta_1 [rad]');
legend('Location', 'Best');
title('System Simulation : \theta_1');
grid on;

subplot(2, 1, 2);
plot(t, rad2deg(q2_response), 'LineWidth', 2, 'DisplayName', '\theta_2');
xlabel('Time [s]');
ylabel('\theta_2 [rad]');
legend('Location', 'Best');
title('System Simulation : \theta_2');
grid on;


% Single Plot
figure;
plot(t, q1_response, 'LineWidth', 2, 'DisplayName', '\theta_1');
hold on;
plot(t, q2_response, 'LineWidth', 2, 'DisplayName', '\theta_2');
xlabel('Time [s]');
ylabel('Angle [rad]');
legend('Location', 'Best');
title('System Simulation : \theta_1 and \theta_2');
grid on;

