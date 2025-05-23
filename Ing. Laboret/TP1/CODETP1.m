% TAREA LABORET 1

close all; clear all; clc;


p1=0;                    %Parametros de G(s)
p2=-3;

k=5;

Tm=0.15;                  %Tiempo de muestreo [s]

G=zpk([],[p1,p2],[k])     %Definicion de G(s)

%FdT DISCRETA de lazo abierto con retentor de orden zero para Tm=0.15s

Gd=c2d(G,Tm,'zoh')
figure(1); pzmap(G); title('Polos y ceros de G(s)'); disp(zero(G)); disp(pole(G));
%discreta:
figure(2); pzmap(Gd); title('Polos y ceros de Gd(z)[Tm:0.15s]'); disp(zero(Gd)); disp(pole(Gd));

%FdT DISCRETA de lazo abierto con retentor de orden zero para Tm=0.9s
Gd1=c2d(G,10*Tm,'zoh')
figure(3);pzmap(Gd1); title('Polos y ceros de Gd(z)[Tm:1.5s]'); disp(zero(Gd1)); disp(pole(Gd1));
%Respuestas al escalon:

figure(4); step(G);
figure(5); step(Gd);
figure(6); step(Gd1);




%Analisis Sistema Discreto:


Kp=dcgain(Gd)
ess=1/(1+Kp)
F=feedback(Gd,1)
figure(7);step(F);grid;

% t=0:Tm:100*Tm % genera rampa
 lsim(F,t,t)
 syms z real
 Kv=((1/Tm)*(z-1)*0.048682*(z-0.8609))/((z-1)*(z-0.6376))
 Kv=limit(Kv,z,1)
 ess=1/Kv

% Generar la entrada tipo rampa
t = 0:Tm:100*Tm; % Vector de tiempo (ajusta el rango según sea necesario)
rampa = t;   % Entrada tipo rampa

% Simular la respuesta del sistema F a la entrada rampa
figure;
lsim(F, rampa, t);
title('Respuesta a una entrada tipo rampa');
xlabel('Tiempo (s)');
ylabel('Salida');
grid on;



%Analisis Lazo cerrado con realimentacion unitaria:

figure(8);rlocus(G);
figure(9);rlocus(Gd);
figure(10);rlocus(Gd1);














