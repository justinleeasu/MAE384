Beta = 0.3;
Lambda = 0.1;
a = 1;
b = 100;
h = 1;
n = (b-a)/h;
T = zeros(1,n+1);
T(1,1) = a;
S = zeros(1,n+1);
I = zeros(1,n+1);
R = zeros(1,n+1);
S(1,1) = 990;
I(1,1) = 10;
R(1,1) = 0;
N = S(1,1)+I(1,1)+R(1,1);

dS = @(S,I) (-Beta/N)*S*I;
dI = @(S,I) (Beta/N)*S*I - Lambda*I;
dR = @(S,I) Lambda*I;

for i = 1:n
    T(i+1) = T(i) + h;
    k1 = dS(S(i),I(i));
    k2 = dS(S(i)+0.5*h,I(i)+0.5*k1*h);
    k3 = dS(S(i)+0.5*h,I(i)+0.5*k2*h);
    k4 = dS(S(i)+h,I(i)+k3*h);
    S(i+1) = S(i) + (1/6)*(k1+2*k2+2*k3+k4)*h;
    k1 = dI(S(i),I(i));
    k2 = dI(S(i)+0.5*h,I(i)+0.5*k1*h);
    k3 = dI(S(i)+0.5*h,I(i)+0.5*k2*h);
    k4 = dI(S(i)+h,I(i)+k3*h);
    I(i+1) = I(i) + (1/6)*(k1+2*k2+2*k3+k4)*h;
    k1 = dR(S(i),I(i));
    k2 = dR(S(i)+0.5*h,I(i)+0.5*k1*h);
    k3 = dR(S(i)+0.5*h,I(i)+0.5*k2*h);
    k4 = dR(S(i)+h,I(i)+k3*h);
    R(i+1) = R(i) + (1/6)*(k1+2*k2+2*k3+k4)*h;
end