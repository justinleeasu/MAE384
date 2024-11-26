Beta = 0.3;
Gamma = 0.1;

a = 1;
b = 100;
h = 1;
n = (b-a)/h;
n = n+1;

T = zeros(1,n+1);
T(1) = a-h;

Susceptible = zeros(1,n+1);
Infected = zeros(1,n+1);
Recovered = zeros(1,n+1);
Susceptible(1) = 990;
Infected(1) = 10;
Recovered(1) = 0;

N = Susceptible(1)+Infected(1)+Recovered(1);

for i = 1:n
    dS = @(T,S) (-Beta/N)*S*Infected(i);
    dI = @(T,I) (Beta/N)*Susceptible(i)*I - Gamma*I;
    dR = @(T,I) Gamma*Infected(i);
    T(i+1) = T(i)+h;
    Susceptible(i+1) = RungeKutta(dS,T(i),Susceptible(i),h);
    Infected(i+1) = RungeKutta(dI,T(i),Infected(i),h);
    Recovered(i+1) = RungeKutta(dR,T(i),Recovered(i),h);
end

function [X] = RungeKutta(dYdX,Xi,Yi,h)
    k1 = dYdX(Xi,Yi);
    k2 = dYdX(Xi+0.5*h,Yi+0.5*k1*h);
    k3 = dYdX(Xi+0.5*h,Yi+0.5*k2*h);
    k4 = dYdX(Xi+h,Yi+k3*h);
    X = Yi + (1/6)*(k1+2*k2+2*k3+k4)*h;
end

figure(1);
plot(T,Susceptible,'r','LineWidth',2);
hold on
plot(T,Infected,'g','LineWidth',2);
plot(T,Recovered,'b','LineWidth',2);
hold off
legend('Susceptible','Infected','Recovered');
