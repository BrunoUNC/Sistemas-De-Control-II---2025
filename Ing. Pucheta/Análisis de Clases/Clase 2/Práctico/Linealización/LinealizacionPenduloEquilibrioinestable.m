clear all; clc; 

syms fi fi_p fi_pp p p_p p_pp M m u long Fricc g; 
disp('Para el equilibrio inestable') 

ang_inic=0; % Punto de equilibrio inestable.

% Aproximaciones válidas para pequeños angulos:
% Cos(fi)~1 ; sin(fi)~fi

p_pp=(1/(M+m))*(u-m*long*fi_pp+m*long*fi_p^2*fi-Fricc*p_p);   %Pequeños angulos 

% fi_pp_=(1/long)*(g*sin(fi)-p_pp*cos(fi));%Expresión completa 
% fi_pp=(1/long)*(g*(fi)-p_pp); %Pequeños angulos para fi~0, sin(fi)~fi, cos(fi)~1 
fi_pp=solve(fi_pp==(1/long)*(g*fi-p_pp),fi_pp); 

%disp('fi_pp='); pretty(simplify(fi_pp)); 
p_pp=subs(p_pp,'fi_pp',fi_pp); 

%disp('p_pp='); pretty(simplify(p_pp)) 


Mat_A=[[0 1 0 0]; 
    [subs(subs(subs(subs(diff(p_pp, p), p,0),p_p,0),fi,ang_inic),fi_p,0), ... 
    subs(subs(subs(subs(diff(p_pp, p_p), p,0),p_p,0),fi,ang_inic),fi_p,0), ... 
    subs(subs(subs(subs(diff(p_pp, fi), p,0),p_p,0),fi,ang_inic),fi_p,0), ... 
    subs(subs(subs(subs(diff(p_pp, fi_p), p,0),p_p,0),fi,ang_inic),fi_p,0)]; 
    [0 0 0 1]; 
    [subs(subs(subs(subs(diff(fi_pp, p), p,0),p_p,0),fi,ang_inic),fi_p,0),... 
    subs(subs(subs(subs(diff(fi_pp, p_p), p,0),p_p,0),fi,ang_inic),fi_p,0),... 
    subs(subs(subs(subs(diff(fi_pp, fi), p,0),p_p,0),fi,ang_inic),fi_p,0),... 
    subs(subs(subs(subs(diff(fi_pp, fi_p), p,0),p_p,0),fi,ang_inic),fi_p,0)]]; 


Mat_B=[0; 
    subs(subs(subs(subs(diff(p_pp, u), p,0),p_p,0),fi,ang_inic),fi_p,0);... 
    0; 
    subs(subs(subs(subs(diff(fi_pp, u), p,0),p_p,0),fi,ang_inic),fi_p,0)]; 

disp('Matriz A:')
pretty(simplify(Mat_A)) 
disp('Matriz B:')
pretty(simplify(Mat_B)) 


