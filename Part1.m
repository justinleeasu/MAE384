Beta = 0.3;
Lambda = 0.1;

a = 1;
b = 100;
h = 1;
n = (b-a)/h;

T = zeros(1,n+1);
T(1) = a;

Susceptible = zeros(1,n+1);
Infected = zeros(1,n+1);
Recovered = zeros(1,n+1);
Susceptible(1) = 990;
Infected(1) = 10;
Recovered(1) = 0;

N = Susceptible(1)+Infected(1)+Recovered(1);


dS = @(S,I) (-Beta/N)*S*I;
dI = @(S,I) (Beta/N)*S*I - Lambda*I;
%dR = @(S,I) Lambda*I;

%Use I(0) to find S(i+1) because I(0) is a constant. So use -B*I(0)/N as a
%constant then dS/dt = C*S(t) and then use Runge-Kutta to find S(i+1). Then
%use dI/dt = C*I(t)-Lambda*I(t) to find B*S(0)/N

for i = 1:n
    T(i+1) = T(i)+h;
    [Susceptible(i+1)] = RungeKutta(dS,Susceptible(i),Infected(i),h);
    [Infected(i+1)] = RungeKutta(dI,Susceptible(i),Infected(i),h);
end

function [X] = RungeKutta(dYdX,Xi,Yi,h)
    k1 = dYdX(Xi,Yi);
    k2 = dYdX(Xi+0.5*k1*h,Yi);
    k3 = dYdX(Xi+0.5*k2*h,Yi);
    k4 = dYdX(Xi+k3*h,Yi);
    X = Xi + (1/6)*(k1+2*k2+2*k3+k4)*h;
end
