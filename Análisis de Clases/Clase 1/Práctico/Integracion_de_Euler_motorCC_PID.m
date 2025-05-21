%Modelado de un motor DC con control PID utilizando método numérico: Integración de Euler

% Desarrollo de método de Euler:
% Inicio: Se comienza con un valor inicial conocido de la variable dependiente (por ejemplo, (y)) en un punto inicial (por ejemplo, (t_0)).
% Tamaño de paso: Se elige un tamaño de paso ((\Delta t)) pequeño. Este es el intervalo de tiempo entre los puntos en los que se estima la solución.
% Cálculo de la pendiente: En cada paso, se calcula la pendiente de la solución en el punto actual. Esta pendiente es dada por la ecuación diferencial que se está resolviendo.
% Estimación del siguiente valor: Se usa la pendiente para estimar el valor de la variable dependiente en el siguiente punto. Específicamente, el nuevo valor ((y_{n+1})) se calcula como (y_n + \Delta t \cdot f(t_n, y_n)), donde (f(t, y)) es la función dada por la EDO.
% Repetición: Este proceso se repite para cada paso sucesivo hasta alcanzar el valor final de (t) deseado.


clear; close all;

X=-[0; 0];     %vector utilizado para guardar las variables de inter�s �ngulo y velocidad
ii=0;          %�ndice de iteraci�n
t_etapa=1e-7;  %Tiempo de muestreo/Paso de integraci�n - Se debe a la din�mica r�pida del sistema      
tF=0.002;      %Tiempo final 
wRef=1;        %Valor de referencia del �ngulo
color_='g';

%Constantes del PID - Seteo
Kp=3; Ki=0.000001; Kd=0.00000007;  %Cuidado, acci�n derivativa causa pico de tensi�n en el inicio
% Kp=1;Ki=0;Kd=0.0001;color_='k';
% Kp=10;Ki=0;Kd=0;color_='b';

Ts=t_etapa; %Tiempo de sampling (muetreo) (Debido a dinámica rápida)

A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts); %C�lculo PID discreto
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;

e=zeros(tF/t_etapa,1);  %Inicializo error como un vector con cantidad de lugares tF/t_etapa en cero
u=0;                    %inicializo acci�n de control

for t=0:t_etapa:tF      % (tF= tiempo final de simulación, se debe a dinámica lenta)
    ii=ii+1; k=ii+2;
    
    X=modmotor(t_etapa, X, u);       % Cálculo de ángulo y velocidad angular a partir del método de Euler 

    e(k)=wRef-X(1);                  % ERROR de la velocidad con respecto a su set point (wRef)
    u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); % Cálculo de la acción de control de un PID
    x1(ii)=X(1);                     % Guardo el valor de Omega (�ngulo) en el vector X para pasarlo a la función modmotor en la próxima iteración
    x2(ii)=X(2);                     % Guardo el valor de wp (velocidad angular) en el vector X para pasarlo a la función modmotor en la próxima iteración
    acc(ii)=u;                       % Guardo el valor de la acción de control
end

t=0:t_etapa:tF;
figure(1);
subplot(3,1,1);hold on;
plot(t,x1,color_);title('Salida y1= �ngulo del motor, \phi');       % Gráfica del �ngulo
subplot(3,1,2);hold on;
plot(t,x2,color_);title('Salida y2= velocidad angular, \Omega_t');  % Gráfica del �ngulo
subplot(3,1,3);hold on;
plot(t,acc,color_);title('Entrada u_t, v_a');  % Gráfica de la acción de control	
xlabel('Tiempo [Seg.]');

%Verificaci�n
 Laa=366e-6;
 J=5e-9;
 Ra=55.6;
 B=0;
 Ki=6.49e-3;
 Km=6.53e-3;
 num=[Ki];
 den=[Laa*J Ra*J+Laa*B Ra*B+Ki*Km ];     %wpp*Laa*J+wp*(Ra*J+Laa*B)+w*(Ra*B+Ki*Km)=Vq*Ki
 sys=tf(num,den)                         %Funci�n de transferencia motor DC velocidad_angular/tensi�n
 figure(2);
 step(sys);title('Respuesta al escal�n');  % Gráfica de la respuesta al escal�n motor solo
 
 


function [X]=modmotor(t_etapa, xant, accion)   %Modela motor DC
    Laa=366e-6; J=5e-9; Ra=55.6; B=0; Ki=6.49e-3; Km=6.53e-3;
    Va=accion;
    h=1e-7; %Permite reducir el paso de integración
    omega= xant(1);
    wp= xant(2);

%Aplicación de la integración de Euler:

    for ii=1:t_etapa/h 
        wpp =(-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+Va*Ki)/(J*Laa); %Derivada de la velocidad angular
        wp=wp+h*wpp; %Cálculo de la velocidad angular a partir de la derivada (aplico Euler) 
        omega = omega + h*wp; %Cálculo del ángulo a partir de la velocidad angular (derivada de la velocidad angular) (aplico Euler)
    end
    X=[omega,wp];  %Guarda en un vector los valores de omega (ángulo) y wp (velocidad angular) calculados por Euler
end
    