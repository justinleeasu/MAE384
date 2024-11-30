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
dS = @(T,S,I) (-Beta/N)*S*I;
dI = @(T,I,S) (Beta/N)*S*I - Gamma*I;
dR = @(T,R,I) Gamma*I;

for i = 1:n
    T(i+1) = T(i)+h;

    k1s = dS(T(i),Susceptible(i),Infected(i));
    k1i = dI(T(i),Infected(i),Susceptible(i));
    k1r = dR(T(i),Recovered(i),Infected(i));

    k2s = dS(T(i)+0.5*h,Susceptible(i)+0.5*k1s*h,Infected(i)+0.5*k1i*h);
    k2i = dI(T(i)+0.5*h,Infected(i)+0.5*k1i*h,Susceptible(i)+0.5*k1s*h);
    k2r = dR(T(i)+0.5*h,Recovered(i)+0.5*k1r*h,Infected(i)+0.5*k1i*h);

    k3s = dS(T(i)+0.5*h,Susceptible(i)+0.5*k2s*h,Infected(i)+0.5*k2i*h);
    k3i = dI(T(i)+0.5*h,Infected(i)+0.5*k2i*h,Susceptible(i)+0.5*k2s*h);
    k3r = dR(T(i)+0.5*h,Recovered(i)+0.5*k2r*h,Infected(i)+0.5*k2i*h);

    k4s = dS(T(i)+h,Susceptible(i)+k3s*h,Infected(i)+k3i*h);
    k4i = dI(T(i)+h,Infected(i)+k3i*h,Susceptible(i)+k3s*h);
    k4r = dR(T(i)+h,Recovered(i)+k3r*h,Infected(i)+k3i*h);

    Susceptible(i+1) = Susceptible(i) + (1/6)*(k1s+2*k2s+2*k3s+k4s)*h;
    Infected(i+1) = Infected(i) + (1/6)*(k1i+2*k2i+2*k3i+k4i)*h;
    Recovered(i+1) = Recovered(i) + (1/6)*(k1r+2*k2r+2*k3r+k4r)*h;
end

figure(1);
plot(T,Susceptible,'r','LineWidth',2);
hold on
plot(T,Infected,'g','LineWidth',2);
plot(T,Recovered,'b','LineWidth',2);
hold off
legend('Susceptible','Infected','Recovered');
