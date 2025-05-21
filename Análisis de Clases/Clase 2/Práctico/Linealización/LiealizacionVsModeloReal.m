clc;clear all; 
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5; % Valores para modelar

%parámetros y valores iniciales
h=0.0001;tiempo=(10/h);p_pp=0;tita_pp=0; t=0:h:tiempo*h; %h:paso de muestreo para integración por método de Euler. 
omega(1)=0;p_p=0:h:tiempo*h; u=linspace(0,0,tiempo+1); 

%Condiciones iniciales. simulación para punto de equilibrio estable
alfa(1)=pi-0-1; color='b';             % Aclaración: alfa es el ángulo del péndulo denominado fi en los scripts anteriores.
%pi-0.1 es una desviación desde el punto de operación?
p(1)=0; p_p(1)=0; u(1)=0; p(1)=0; i=1;

%Versión linealizada en el equilibrio estable. Sontag Pp 104. 
% Vector de estado=[p(i); p_p(i); alfa(i); omega(i)]   cuidado, omega es fi_p en los scripts anteriores.
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(long*M) -g*(m+M)/(long*M) 0] 
Mat_B=[0; 1/M; 0; 1/(long*M)]                                                    
X0=[0 0 pi 0]';x=[0 0 alfa(1) 0]'; 

while(i<(tiempo+1)) 
    
    %Variables del sistema no lineal 
    estado=[p(i); p_p(i); alfa(i); omega(i)]; 
    u(i)=0; 
    %Sistema no lineal (serían f1, f2, f3 y f4)
    p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))-Fricc*p_p(i)); 
    tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i))); 
    p_p(i+1)=p_p(i)+h*p_pp; 
    p(i+1)=p(i)+h*p_p(i); 
    omega(i+1)=omega(i)+h*tita_pp; 
    alfa(i+1)=alfa(i)+h*omega(i); 
    
    
    %Variables del sistema lineal 
    pl(i)=x(1); p_pl(i)=x(2);alfal(i)=x(3);omegal(i)=x(4); 
    %Sistema lineal (Linealización a partir de f1, f2, f3, f4)
    xp=Mat_A*(x-X0)+Mat_B*u(i); 
    x=x+h*xp; 
    i=i+1; 
    
end 


pl(i)=x(1); p_pl(i)=x(2);alfal(i)=x(3);omegal(i)=x(4); 
figure(1);hold on; 
subplot(3,2,1);plot(t,omega,color);grid on; title('Velocidad Ángulo');
hold on;plot(t,omegal,'k'); 
subplot(3,2,2);plot(t,alfa,color);hold on;
plot(t,pi*ones(size(t)),'k');plot(t,alfal,'k');
grid on;title('Ángulo');hold on;
subplot(3,2,3); plot(t,p,color);grid on;title('Posición carro');
hold on;plot(t,pl,'k');
subplot(3,2,4);plot(t,p_p,color);grid on;title('Velocidad carro');
hold on;plot(t,p_pl,'k');
subplot(3,1,3);plot(t,u,color);grid on;title('Acción de control');xlabel('Tiempo en Seg.');
hold on; 
figure(2);hold on;
subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Ángulo');ylabel('Velocidad angular');hold on;
subplot(2,2,1);plot(alfal,omegal,'k');
subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posición carro');ylabel('Velocidad carro');hold on;
subplot(2,2,2);plot(pl,p_pl,'k');





