%% Caso de estudio 3. Sistema no lineal de cuatro variables de estado con y sin observador (todo en tiempo continuo)
%% Desplazamiento hasta 10m

clear all;
close all;
clc;

%% Parámetros físicos del sistema
m = 0.1;    % Masa del péndulo [kg]
F = 0.1;    % Coeficiente de fricción [N·s/m]
l = 1.6;    % Longitud del péndulo [m]
g = 9.8;    % Aceleración gravitatoria [m/s²]
M = 1.5;    % Masa del carro [kg]

%% Matrices del sistema continuo linealizado alrededor del equilibrio (delta = 0, phi = 0)
A = [0 1 0 0;
     0 -F/M -m*g/M 0;
     0 0 0 1;
     0 -F/(l*M) -g*(m+M)/(l*M) 0];
B = [0; 1/M; 0; 1/(l*M)];
C = [1 0 0 0; 0 0 1 0];  % Salidas: desplazamiento y ángulo

%% Configuración temporal
T = 25;      % Tiempo total de simulación [s]
h = 1e-4;    % Paso de integración [s]
t = 0:h:T;

%% Ampliación del sistema para control integral (solo en delta)
Cref = C(1, :);
aA = [A, zeros(4, 1); -Cref*A, 0];
aB = [B; -Cref*B];
Q = diag([10 0.1 10 0.1 3]);
R = 2;
K = lqr(aA, aB, Q, R);
Kp = K(1:4);
Ki = -K(5);

%% Observador de Luenberger en tiempo continuo (estimación de todo el estado a partir de delta y phi)
Ao = A';
Bo = C';
Qo = diag([10 1 10 1]);
Ro = diag([10 10]);  % <- Matriz R ahora es diagonal (dual al número de salidas)
Ko = lqr(Ao, Bo, Qo, Ro)';

%% Inicialización de estados y señales
n = length(t);
delta = zeros(1, n); deltaP = zeros(1, n); phi = zeros(1, n); omega = zeros(1, n);
delta_hat = zeros(1, n); deltaP_hat = zeros(1, n); phi_hat = zeros(1, n); omega_hat = zeros(1, n);
u_real = zeros(1, n); u_obs = zeros(1, n);
reference = 10 * ones(1, n);

x_real = [0; 0; pi; 0];
x_hat = [0; 0; 0; 0];
x_op = [0; 0; pi; 0];
e_int = 0;
e_int_hat = 0;

%% Simulación en tiempo continuo
for k = 1:n-1
    y = C * x_real;
    y_hat = C * (x_hat - x_op);

    err = reference(k) - y(1);
    err_hat = reference(k) - y_hat(1);
    e_int = e_int + h * err;
    e_int_hat = e_int_hat + h * err_hat;

    % Control LQR con acción integral
    u_k = -Kp * (x_real - x_op) + Ki * e_int;
    u_hat_k = -Kp * (x_hat - x_op) + Ki * e_int_hat;

    % Zona muerta
   % if abs(u_k) < 1, u_k = 0; else, u_k = sign(u_k) * (abs(u_k) - 1); end
   % if abs(u_hat_k) < 1, u_hat_k = 0; else, u_hat_k = sign(u_hat_k) * (abs(u_hat_k) - 1); end

    % Dinámica real no lineal
    phi_pp = (1/l)*(g*sin(x_real(3)) - ((1/(M+m))*(u_k - m*l*0*cos(x_real(3)) + m*l*x_real(4)^2*sin(x_real(3)) - F*x_real(2))) * cos(x_real(3)));
    x_real(2) = x_real(2) + h*(1/(M+m))*(u_k - m*l*phi_pp*cos(x_real(3)) + m*l*x_real(4)^2*sin(x_real(3)) - F*x_real(2));
    x_real(1) = x_real(1) + h * x_real(2);
    x_real(4) = x_real(4) + h * phi_pp;
    x_real(3) = x_real(3) + h * x_real(4);

    % Observador (estimación del estado)
    x_hat_dot = A*(x_hat-x_op) + B*u_hat_k + Ko*(y - y_hat);
    x_hat = x_hat + h * x_hat_dot;

    % Guardar señales
    delta(k+1) = x_real(1); deltaP(k+1) = x_real(2); phi(k+1) = x_real(3); omega(k+1) = x_real(4);
    delta_hat(k+1) = x_hat(1); deltaP_hat(k+1) = x_hat(2); phi_hat(k+1) = x_hat(3); omega_hat(k+1) = x_hat(4);
    u_real(k) = u_k; u_obs(k) = u_hat_k;
end

%% Gráficas comparativas
figure;
subplot(3,1,1);
plot(t, delta, 'b', t, delta_hat, 'r--'); grid on;
title('Desplazamiento \delta'); ylabel('[m]'); legend('Real','Estimado');

subplot(3,1,2);
plot(t, phi, 'b', t, phi_hat, 'r--'); grid on;
title('Ángulo \phi'); ylabel('[rad]'); legend('Real','Estimado');

subplot(3,1,3);
plot(t, u_real, 'b', t, u_obs, 'r--'); grid on;
title('Acción de control'); ylabel('[V]'); xlabel('Tiempo [s]'); legend('Sin observador','Con observador');
