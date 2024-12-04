b = [30 10];

for j = 1:length(b)
    clear y
    Beta = 0.3;
    Gamma = 0.1;

    a = 1;
    h = 1;
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
    dS = @(T,S,I) (-Beta/N)*S*I;
    dI = @(T,I,S) (Beta/N)*S*I - Gamma*I;
    dR = @(T,R,I) Gamma*I;
    y(1) = log(Infected(1));

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
        y(i+1) = log(Infected(i+1));
    end

    A1NumPart1 = length(T)*sum(T.*y);
    A1NumPart2 = (sum(T))*(sum(y));
    A1Num = A1NumPart1-A1NumPart2;

    A1DenomPart1 = length(T)*sum((T.^2));
    A1DenomPart2 = (sum(T))^2;
    A1Denom = A1DenomPart1-A1DenomPart2;
    
    A1 = A1Num/A1Denom;
    A0Part1 = (1/length(T))*sum(y);
    A0Part2 = A1*(1/length(T))*sum(T);
    A0 = A0Part1 - A0Part2;

    I0(j) = exp(A0);
    B(j) = (A1+Gamma)*(N/Susceptible(1));
end

fprintf('The estimation of initial infected for T = 30 is %4.4f\n', I0(1));
fprintf('The estimation of initial infected for T = 10 is %4.4f\n', I0(2));
fprintf('The estimation of beta for T = 10 is %4.4f\n', B(1));
fprintf('The estimation of beta for T = 10 is %4.4f\n', B(2));
fprintf('The true value for initial infected is %4.1f\n',Infected(1));
fprintf('The true value for beta is %4.1f\n',Beta);
fprintf('\nClearly, the estimates get better with a smaller set of data.\nThis is because the smaller step size means that the linear\nleast squares regression model will be more accurate as it has\nless data to average together.');
