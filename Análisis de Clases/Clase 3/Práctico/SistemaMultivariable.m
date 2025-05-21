%VerificaciÃ³n de la EDO 

%Sistema de ecuaciones diferenciales MIMO a trabajar:
% y1pp+y1p+y2=u1+u2p
% y2p=y1+u1p,

% Condic iniciales nulas, respuestas al escalon.

clear all;clc;close all;

T=10; Kmax=20000; At=T/Kmax; TamanioFuente=12;
u1=zeros(1,Kmax); y1=u1; y2=u1; y1p=y1; y2p=y2;
t=0:At:T-At;y1pp=0;

u1=[0 2*cos(10*t)]; u1p=[diff(u1)/At];   % Señales de entrada y sus derivadas
u2=[0 2*sin(1*t)]; u2p=[diff(u2)/At];

for jj=1:Kmax-1                          % EDO + método de integración por Euler
% y1pp+y1p+y2=u1+u2p;
% y2p=y1+u1p;
 y1pp=-y1p(jj)-y2(jj)+u1(jj)+u2p(jj); 
 y2p(jj)=y1(jj)+u1p(jj);
 %Integro
 y2(jj+1) =y2(jj) +y2p(jj)*At;
 y1p(jj+1)=y1p(jj)+y1pp*At;
 y1(jj+1) =y1(jj) +y1p(jj)*At;
end

figure;title('Sistema MIMO');             %Plots EDOs
subplot(2,2,1);plot(t,y1,'.k');title('y_1'),hold on        
subplot(2,2,2);plot(t,y2,'.k');title('y_2'); hold on
subplot(2,2,3);plot(t,u1(1:numel(t)),'k');title('u_1');xlabel('t [seg.]');
subplot(2,2,4);plot(t,u2(1:numel(t)),'k');title('u_2');xlabel('t [seg.]');



%Modelado en variables de estado:
A=[0,1,0 ; 0,-1,-1 ; 1,0,0];
B=[0, 1 ; 0, -1 ; 0,0];
C=[1, 0,0 ; 0,0,1];
D=[0 0 ; 1 0];
Sis_ve=ss(A,B,C,D)                        %sistema matricial en variables de estado
[num, den]=ss2tf(Sis_ve.a,Sis_ve.b,Sis_ve.c,Sis_ve.d,1) %Cálculo de funciones de transferencia con respecto a entrada 1 (y1/u1 y2/u1)
%sys1=tf(num,den) 
[num, den]=ss2tf(Sis_ve.a,Sis_ve.b,Sis_ve.c,Sis_ve.d,2) %Cálculo de funciones de transferencia con respecto a entrada 1 (y1/u2 y2/u2)
%sys2=tf(num,den)

x=[0;0;0];U=[0;0];y_1(1)=0;y_2(1)=0; %Valores iniciales
for jj=1:Kmax-1                           %Sistema en variables de estado + integración por Euler
 U(1)=u1(jj);
 U(2)=u2(jj);
 xp=A*x+B*U;
 y =C*x+D*U;
 x=x+xp*At;
 y_1(jj+1)=y(1);
 y_2(jj+1)=y(2);
end
%Plots sistema en variables de estado:
subplot(2,2,1);plot(t,y_1,'r');legend('EDO','Evolución en Varaibles de Estado');         
subplot(2,2,2);plot(t,y_2,'r');legend('EDO','Evolución en Varaibles de Estado');
legend('boxoff');title('y_2','FontSize',TamanioFuente);