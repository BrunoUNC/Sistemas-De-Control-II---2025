%VerificaciÃ³n de la EDO Ejemplo

% Caso de ejemplo:  ypp+yp+y=u1+up

% Condiciones iniciales nulas, respuestas al seguimiento.
clear all;close all;clc;

T=10; Kmax=10000; At=T/Kmax; TamanioFuente=12;
u=zeros(1,Kmax); y=u;yp=y; t=0:At:T-At;
u=[0 sin(5*t+2)]; up=[0 diff(u)/At];      % Entrada senoidal



for jj=1:Kmax-1                           % Método de integración de Euler aplicado
 % ypp+yp+y=u1+u1p                          al sistema real (EDO)
 ypp=-yp(jj)-y(jj)+u(jj)+up(jj);
 %Integro
 yp(jj+1)=yp(jj)+ypp*At;
 y(jj+1) =y(jj) +yp(jj)*At;
end

figure (1);title('Sistema SISO');
subplot(2,2,1);plot(t,y,'.k');title('y'),hold on                             %plot de la salida de la EDO
subplot(2,2,3);plot(t,u(1:numel(t)),'k');title('u');xlabel('timepo [seg.]'); %plot de la entrada


%Modelado en variables de estado:
A=[0,1;-1,-1];
B=[1;0];
C=[1, 0];
D=[0];
Sis_ve=ss(A,B,C,D) %sistema en variables de estado
[num, den]=ss2tf(Sis_ve.a,Sis_ve.b,Sis_ve.c,Sis_ve.d,1) %cálculo de función de transferencia
x=[0;0];
U=[0]; y_1(1)=0; hh=At*1; tve(1)=0;

for jj=1:Kmax-1                         % Método de integración de Euler aplicado
 U=u(jj);                               % al sistema en variables de estado
 xp=A*x+B*U;
 Y =C*x+D*U;
 x=x+xp*At;
 y_1(jj+1)=Y(1);
 tve(jj+1)=tve(jj)+hh;                  % vector de tiempo
end

subplot(2,2,1);plot(tve,y_1,'r');                                            %plot de la salida en variable de estado
legend('EDO','Evolución en Variables de Estado');                                      
legend("boxoff");