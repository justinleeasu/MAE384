b = [30 10];

for j = 1:length(b)

    Beta = 0.3;
    Gamma = 0.1;

    h = 1;
    a = 1;
    n = (b(j)-a)/h;
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
    y = zeros(1,n+1);
    y(1,1) = log(Infected(1));

    for i = 1:n
        dS = @(T,S) (-Beta/N)*S*Infected(i);
        dI = @(T,I) (Beta/N)*Susceptible(i)*I - Gamma*I;
        dR = @(T,I) Gamma*Infected(i);
        T(i+1) = T(i)+h;
        Susceptible(i+1) = RungeKutta(dS,T(i),Susceptible(i),h);
        Infected(i+1) = RungeKutta(dI,T(i),Infected(i),h);
        Recovered(i+1) = RungeKutta(dR,T(i),Recovered(i),h);
        y(i+1) = log(Infected(i+1));
    end

    A1Num = n*sum((T.*y)) - sum(T)*sum(y);
    A1Denom = n*sum((T.^2)) - (sum(T))^2;
    A1 = A1Num/A1Denom;

    A0 = (1/(n))*sum(y) - A1*(1/(n))*sum(T);
    I0(j) = exp(A0);
    B(j) = (A1+Gamma)*(N/Susceptible(1));
    
    clear x
    syms x
end

function [X] = RungeKutta(dYdX,Xi,Yi,h)
k1 = dYdX(Xi,Yi);
k2 = dYdX(Xi+0.5*h,Yi+0.5*k1*h);
k3 = dYdX(Xi+0.5*h,Yi+0.5*k2*h);
k4 = dYdX(Xi+h,Yi+k3*h);
X = Yi + (1/6)*(k1+2*k2+2*k3+k4)*h;
end
