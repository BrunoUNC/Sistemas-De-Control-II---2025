% Bruno Gonzalez - Sistemas No Lineales

clear all; clc;

%% Parámetros físicos del sistema
m = 1;              % masa [kg]
b = 0.4;            % fricción [N.m.s]
l = 1;              % longitud [m]
g = 10;             % gravedad [m/s^2]
delta_deg = 90;     % ángulo de referencia en grados
delta = deg2rad(delta_deg);  % en radianes

% Torque de equilibrio
uf = m * g * l * sin(delta);  % T_eq = mgl * sin(delta)

%% Linealización del modelo en el punto de equilibrio
try
    [A, B, C, D] = linmod('pendulo_mod_tarea', [0; 0], uf);  %vector X=[0;0] ya que se x1 desplaza el equilibrio al origen
catch
    error('No se pudo linealizar el modelo. Verifica que "pendulo_mod_tarea.slx" esté correctamente configurado.');
end
% Asegurate de tener el modelo 'pendulo_mod_tarea.slx' abierto o accesible

%% Resultados
disp('--- Matriz A ---');
disp(A);

disp('--- Matriz B ---');
disp(B);

disp('--- Matriz C ---');
disp(C);

disp('--- Matriz D ---');
disp(D);

%% Análisis de estabilidad (Lyapunov indirecto)
autovalores = eig(A);
disp('--- Autovalores del sistema linealizado ---');
disp(autovalores);

if all(real(autovalores) < 0)
    disp('El sistema linealizado es localmente asintóticamente estable.');
elseif any(real(autovalores) > 0)
    disp('El sistema linealizado es inestable.');
else
    disp('El método indirecto no permite concluir la estabilidad (autovalores en j?).');
end

%% Verificación de controlabilidad
Co = ctrb(A, B);
rango = rank(Co);
disp(['--- Rango de la matriz de controlabilidad: ', num2str(rango), ' de ', num2str(size(A,1))]);

if rango == size(A,1)
    disp('El sistema es completamente controlable.');
else
    disp('El sistema NO es completamente controlable.');
end


%% Diseño de Controlador con Acción Integral para el Péndulo Linealizado

%% Definición del sistema
A = [0 1; -10 -0.4];
B = [0; 1];
C = [1 0];

%% Construcción del sistema extendido Aa y Ba
Aa = [[A; C], zeros(3,1)];  % A extendido: A con fila C y columna de ceros
Ba = [B; 0];                % B extendido: agregando un 0

%% Verificación de controlabilidad y estabilidad
disp('Autovalores del sistema extendido Aa:');
disp(eig(Aa));

disp('Rango de la matriz de controlabilidad extendida:');
disp(rank(ctrb(Aa, Ba)));

%% Diseño del controlador por asignación de polos
p = -2;                            % Polo deseado triple
K = acker(Aa, Ba, [p p p]);        % Cálculo de la ganancia 

k1 = K(1);
k2 = K(2);
k3 = K(3);

disp('Ganancia del controlador:');
disp(['k1 = ', num2str(k1)]);
disp(['k2 = ', num2str(k2)]);
disp(['k3 = ', num2str(k3)]);

%% Verificación de la dinámica del sistema en lazo cerrado
disp('Autovalores del sistema en lazo cerrado:');
disp(eig(Aa - Ba*K));

%% Estimación del tiempo de respuesta (regla práctica)
tscalc = 7.5 / (-p);
disp(['Tiempo de respuesta estimado: ', num2str(tscalc), ' segundos']);

%% Simulación del sistema controlado con acción integral
sim('pendulo_pid_tarea');

%% Gráfica 1: Salida del sistema
figure(1);
plot(tout, yout, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Tiempo [s]');
ylabel('Ángulo del péndulo');
title('Salida del sistema');

%% Gráfica 2: Plano de fases (x1 vs x2)
figure(2);
plot(yout, velocidad, 'r', 'LineWidth', 1.5);
grid on;
xlabel('x_1(t) (posición)');
ylabel('x_2(t) (velocidad)');
title('Plano de fases');

%% Gráfica 3: Torque aplicado
figure(3);
plot(tout, torque, 'k', 'LineWidth', 1.5);
grid on;
xlabel('Tiempo [s]');
ylabel('Torque u(t)');
title('Torque total aplicado');

%% Gráfica 4: Acción integral
figure(4);
plot(tout, -accint, 'm', 'LineWidth', 1.5);
grid on;
xlabel('Tiempo [s]');
ylabel('Acción integral');
title('Acción integral acumulada');

%% Métricas de desempeño
delta = 90; % referencia final deseada

% Sobrepaso máximo
ymax = max(yout);
S = (ymax - delta) / delta * 100;

% Error relativo instantáneo
erel = (delta - yout) / delta;
efinal = erel(end);  % error final en régimen permanente

% Tiempo de establecimiento: último instante con error > 2%
ind = find(abs(erel) > 0.02);
if ~isempty(ind)
    tss = tout(ind(end));      % tiempo de establecimiento
    yte = yout(ind(end));      % salida en ese instante
else
    tss = NaN;
    yte = NaN;
end

% Valores finales de control
uf = torque(end);              % torque final
Intf = -accint(end);           % acción integral final

%% Mostrar resultados
fprintf('--- Resultados de la simulación ---\n');
fprintf('Sobrepaso: %.2f %%\n', S);
fprintf('Error final relativo: %.4f\n', efinal);
fprintf('Tiempo de establecimiento (±2%%): %.2f s\n', tss);
fprintf('Salida al tiempo de establecimiento: %.4f °\n', yte);
fprintf('Torque final aplicado: %.4f N·m\n', uf);
fprintf('Acción integral final: %.4f\n', Intf);

%Recordar variar m en simulink para verificar robustez