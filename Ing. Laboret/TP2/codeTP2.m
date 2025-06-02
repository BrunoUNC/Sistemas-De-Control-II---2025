clc; clear; close all;

%% 1. Datos del sistema
G_s = tf(5, [1 3 0]);    % G(s) = 5 / (s*(s+3))
T = 0.15;               % Tiempo de muestreo

%% 2. Discretización con ZOH
Gd = c2d(G_s, T, 'zoh');
disp('Planta discretizada Gd(z):');
Gd

%% 3. Especificaciones de diseño
sobrepaso = 5;         % %
ts = 2;                % tiempo de establecimiento (2%)
xi = -log(sobrepaso/100)/sqrt(pi^2 + log(sobrepaso/100)^2);   % amortiguamiento
wn = 4/ts;             % frecuencia natural
wd = wn*sqrt(1 - xi^2); % frecuencia amortiguada

%% 4. Polos deseados en el plano z
r = exp(-xi*wn*T);         % módulo
theta = wd*T;              % argumento
z1 = r*exp(1j*theta);      % polo 1
z2 = r*exp(-1j*theta);     % polo 2

fprintf('\nPolos deseados en z:\n z1 = %.4f + %.4fj\n z2 = %.4f %.4fj\n', real(z1), imag(z1), real(z2), imag(z2));

% Número de muestras por ciclo
m = 2*pi / (wd*T);
fprintf('Número de muestras por ciclo: m = %.2f\n', m);

%% 5. Lugar de raíces y diseño del controlador
disp('Abrir Sisotool para diseñar el controlador:');
sisotool(Gd);

% Una vez diseñado el controlador C, exportar al workspace desde sisotool

% Supongamos que se exportó como C
% (descomentar las siguientes líneas luego de exportar desde sisotool)

 F = feedback(C*Gd,1);
% 
% %% 6. Verificación de la respuesta
 figure;
 step(F);
 title('Respuesta al escalón del sistema en lazo cerrado');
% 
% % Error en régimen permanente
 ess = 1 - dcgain(F);
 fprintf('\nError en régimen permanente ante escalón: %.5f\n', ess);
% 


fprintf('\nFunción de transferencia a lazo cerrado F(z):\n');
F

fprintf('\Compensador C):\n');
C

%% Mostrar polos y ceros
fprintf('\nPolos del sistema cerrado:\n');
disp(pole(F));

fprintf('\nCeros del sistema cerrado:\n');
disp(zero(F));

% % Mapa de polos y ceros
 figure;
 pzmap(F);
 title('Polos y ceros del sistema a lazo cerrado');
 
 %% 7. Cálculo de desempeño dinámico
[y, t] = step(F);
y_final = y(end);
y_max = max(y);
overshoot = (y_max - y_final) / y_final * 100;

idx_ts = find(abs(y - y_final) > 0.02 * y_final, 1, 'last');
if isempty(idx_ts)
    ts2 = 0;
else
    ts2 = t(idx_ts);
end

fprintf('\nSobrepaso: %.2f %%\n', overshoot);
fprintf('Tiempo de establecimiento (2%%): %.3f s\n', ts2);