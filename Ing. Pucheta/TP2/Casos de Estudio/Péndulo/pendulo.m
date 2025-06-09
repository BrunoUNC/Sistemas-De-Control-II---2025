%% Caso de estudio 3. Sistema no lineal de cuatro variables de estado. Péndulo Grúa.
%% Sin observador

clear all;
close all;
clc;

%% Definición de parámetros del sistema
m = 0.1;    % Masa del péndulo [kg]
F = 0.1;    % Coeficiente de fricción [N.s/m]
l = 1.6;    % Longitud del péndulo [m]
g = 9.8;    % Aceleración gravitatoria [m/s^2]
M = 1.5;    % Masa del carro [kg]

%% Matrices del sistema continuo linealizado alrededor de (delta=0, phi=0)
Ac = [0 1 0 0;
      0 -F/M -m*g/M 0;
      0 0 0 1;
      0 -F/(l*M) -g*(m+M)/(l*M) 0];
Bc = [0; 1/M; 0; 1/(l*M)];
Cc = [1 0 0 0; 0 0 1 0];  % Salidas: desplazamiento y ángulo
Dc = [0];

%% Tiempos de muestreo y simulación
ts = 1e-2;      % Tiempo de muestreo [s]
T = 25;        % Tiempo total de simulación [s]
h = 1e-4;      % Paso de integración [s]
Kmax = T / ts; % Cantidad de iteraciones del lazo externo

%% Discretización del sistema
sys = ss(Ac, Bc, Cc, Dc);
dSys = c2d(sys, ts, 'zoh');
A = dSys.A;
B = dSys.B;

%% Ampliación del sistema para referencia en desplazamiento (delta)
Cref = Cc(1,:);  % Solo se considera la salida del desplazamiento
aA = [A, zeros(4,1); -Cref*A, 1];
aB = [B; -Cref*B];
Q = diag([0.1 1e-2 1 0.1 9.3e-6]);
R = 1.9e-4;
K = dlqr(aA, aB, Q, R);
Kp = K(1:4);     % Parte proporcional (retroalimentación de estados)
Ki = -K(5);      % Parte integral (integrador del error)

%% Inicialización de variables
phi(1) = pi;                      % Condición inicial del ángulo
x = [0; 0; phi(1); 0];            % Estado inicial del sistema
x_op = [0; 0; pi; 0];             % Punto de operación deseado
v(1) = 0;                         % Estado del integrador del error
reference(1) = 10;                % Referencia deseada para delta

%% Simulación del sistema no lineal
K = Kp;
KI = Ki;
i = 1;

delta = x(1);
deltaP = x(2);
omega = x(4);

for index = 1:Kmax
    y = Cc*x;
    v(index+1) = v(index) + reference(index) - y(1);
    u1(index) = -K*(x - x_op) + KI*v(index+1);

    % Zona muerta en la acción de control
    if abs(u1(index)) < 1
        u1(index) = 0;
    else
        u1(index) = sign(u1(index)) * (abs(u1(index)) - 1);
    end

    for j = 1:(ts/h)
        % Dinámica no lineal del sistema
        u(i) = u1(index);
        phiPP = (1/l)*(g*sin(x(3)) - ((1/(M+m))*(u(i) - m*l*0*cos(x(3)) + m*l*x(4)^2*sin(x(3)) - F*x(2))) * cos(x(3)));
        x(2) = x(2) + h*(1/(M+m))*(u(i) - m*l*phiPP*cos(x(3)) + m*l*x(4)^2*sin(x(3)) - F*x(2));
        x(1) = x(1) + h*x(2);
        x(4) = x(4) + h*phiPP;
        x(3) = x(3) + h*x(4);

        % Guardado para gráficas
        delta(i+1) = x(1);
        deltaP(i+1) = x(2);
        phi(i+1) = x(3);
        omega(i+1) = x(4);

        reference(i+1) = reference(i);
        i = i+1;
    end
end

u(i) = u1(index);
t = 0:h:(length(delta)-1)*h;

%% Gráficas finales en una única ventana
figure;
subplot(5,1,1);
plot(t, delta, 'LineWidth', 1.5); hold on;
plot(t, reference(1:length(t)), 'r--', 'LineWidth', 1.5); hold off; grid on;
xlabel('Tiempo [seg]'); ylabel('Distancia [m]'); title('Desplazamiento \delta'); legend('\delta', 'Referencia');

subplot(5,1,2);
plot(t, deltaP, 'LineWidth', 1.5); grid on;
xlabel('Tiempo [seg]'); ylabel('Velocidad [m/s]'); title('Velocidad de desplazamiento');

subplot(5,1,3);
plot(t, phi, 'LineWidth', 1.5); grid on;
xlabel('Tiempo [seg]'); ylabel('Ángulo [rad]'); title('Ángulo \phi');

subplot(5,1,4);
plot(t, phi/pi*100 - 100, 'LineWidth', 1.5); grid on;
xlabel('Tiempo [seg]'); ylabel('Error [%]'); title('Ángulo \phi respecto a \pi');

subplot(5,1,5);
plot(t, omega, 'LineWidth', 1.5); grid on;
xlabel('Tiempo [seg]'); ylabel('Velocidad angular [rad/s]'); title('Velocidad angular');

figure;
plot(phi, omega, 'LineWidth', 1.5); grid on;
xlabel('Ángulo \phi [rad]'); ylabel('Velocidad angular \omega [rad/s]'); title('Plano de fase: \phi vs \omega');

figure;
plot(delta, deltaP, 'LineWidth', 1.5); grid on;
xlabel('Distancia \delta [m]'); ylabel('Velocidad \delta'' [m/s]'); title('Plano de fase: \delta vs \delta''');
