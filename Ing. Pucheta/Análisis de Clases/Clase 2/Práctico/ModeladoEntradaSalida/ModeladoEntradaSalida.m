%Modelado a partir de resupuesta al escal�n.
%Respuesta de segundo orden con polos complejos
%cojugados

clc; clear all;
% Kuo Sistemas de control autom�tico Fig 7-18, Secci�n 7-5 Pp395.
%Codigo realizado por JAP 

wn_r=1.5;sita_r=.85;TMax=20;
% wn_r=1.5;sita_r=.5;TMax=50;
% wn_r=1.5;sita_r=.15;TMax=120;


% FdT real:
sys_G=tf([wn_r^2],[1 2*sita_r*wn_r wn_r^2]);
sys_real_=sys_G

%Respuesta al escal�n:
[y,t0]=step(sys_G,TMax/wn_r);
step(sys_G,TMax/wn_r)

%Pubntos de muestras: m�ximo y m�nimo.
[ymax a ]=max(y);
t1=t0(a);
[ymin b ]=min(y(a:end));
t2=t0(b+a-1);
hold on

plot(t1,ymax,'pk');
plot(t2,ymin,'sk');

%C�luclo de wn a partir de la frecuencia detectada por los puntos m�ximos y
%m�nimos de la din�mica. wn=2pi.f=2pi/T=2pi/(2(t2-t1))
w_n=(2*pi)/(2*(t2-t1));

%C�lculo de sistema identificado:
beta=-log(ymax-1)/pi;

sita_id=beta/(sqrt(1+beta^2));

wn_id=w_n/(sqrt(1-sita_id^2));

den_id=[1 2*sita_id*wn_id wn_id^2];

sys_id=tf(wn_id^2,den_id)
step(sys_id)
legend('Real','M�ximo','Valle','Identificada')
