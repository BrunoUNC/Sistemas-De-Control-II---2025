clc; clear;

%% 1. Planta continua y discretizaci�n
G = tf(5, [1 3 0]);
T = 0.15;
Gd = c2d(G, T, 'zoh');

fprintf('Planta discretizada Gd(z):\n');
Gd

%% 2. Abrir sisotool para dise�ar el controlador PI
% Desde esta interfaz se agregar� un cero y un polo en z=1
sisotool(Gd)

% Una vez dise�ado el PI y exportado como 'C', continuar:
% (Descomentar luego de exportar)

 F = feedback(C*Gd, 1);

 fprintf('\nFunci�n de transferencia del sistema cerrado con PID:\n');
 F
 
 %% Respuesta al escal�n y an�lisis
[y, t] = step(F);
y_final = y(end);

% Sobrepaso
y_max = max(y);
overshoot = (y_max - y_final) / y_final * 100;

% Tiempo de establecimiento al 2%
% Definido como el tiempo a partir del cual la salida se mantiene dentro del �2% de y_final
idx_settling = find(abs(y - y_final) > 0.02 * y_final, 1, 'last');
if isempty(idx_settling)
    ts_2 = 0;
else
    ts_2 = t(idx_settling);
end

%% Resultados
fprintf('Sobrepaso: %.2f %%\n', overshoot);
fprintf('Tiempo de establecimiento (2%%): %.3f segundos\n', ts_2);

 fprintf('\nPolos del sistema cerrado:\n'); disp(pole(F));

 fprintf('\nCeros del sistema cerrado:\n');
 disp(zero(F));

 figure;
 pzmap(F);
 title('Polos y ceros del sistema con PID');

 figure;
 step(F);
 title('Respuesta al escal�n del sistema con PID');

 ess = 1 - dcgain(F);
 fprintf('\nError en r�gimen permanente ante escal�n: %.5f\n', ess);