function asa_da_vinci(w1)

r1 = 1.35;
r2 = 3.9;

w2 = (-w1*r1)/r2; %Vel angular engrenagem maior


%Tamanho das barras do mecanismo

R = sqrt(2);
L1 = sqrt(20);
L2 = sqrt(66);
L3 = 8.3;
L4 = 18.1;
L5 = 29.2;
L6 = 2.7;
d = 2;

%Números de pontos a ser analizados dentro do espectro da lista que varia
%de 0 a 360.
numPontos = 100;
Theta = linspace(0,360,numPontos);

%Matrizes que armazenam os valores das velocidades e posições
X1 = zeros(numPontos,1);
X2 = zeros(numPontos,1);
VP3 = zeros(numPontos,2);
VelP4 = zeros(numPontos,2);
VP6 = zeros(numPontos,2);
VP5 = zeros(numPontos,2);
VP7 = zeros(numPontos,2);
VP8 = zeros(numPontos,2);
VP9 = zeros(numPontos,2);

%Matrizes que armazenam os valores das velocidades angulares
W3 = zeros(numPontos,1);
W4 = zeros(numPontos,1);
W5 = zeros(numPontos,1);
W6 = zeros(numPontos,1);
W7 = zeros(numPontos,1);
W8 = zeros(numPontos,1);


%Definições para ode45
f = abs(w2/2/pi);
tspan = [ 0 1/f ];
rB0 = [ r2 0 ];

%Cálculo das posições do ponto P3
[T, PosP3] = ode45(@(t,PosP3)odeFun(t,PosP3,w2,r2),tspan,rB0);


for i = 1: numPontos
    theta = Theta(i)*pi/180; %Transformando o theta em radianos
    
    senot = sin(theta);
    cost = cos(theta);

    %Velocidades no ponto P3
    VP3(i,1) = -w2*r2*senot; 
    VP3(i,2) = w2*r2*cost;
    
    %Cálculo do sistema não linear
    if i == 1
        X0 = [0 0];
    else
        X0 = [X1(i-1) X2(i-1)];
    end
    [x, ~] = fsolve(@(x)funcaoMecanismoDeLee(x,theta,L1,L2,R),X0);
    %Resultados do sistema não linear (sinGamma e sinAlpha)
    X1(i) = x(1);
    X2(i) = x(2);
    sinGamma = x(1);
    sinAlpha = x(2);
    cosGamma = sqrt(1-sinGamma^2);
    cosAlpha = sqrt(1-sinAlpha^2);
    
    %Cálculo das velocidades no ponto P4
    PosP4(i,1) = R*cos(theta) - L1*cosGamma;
    PosP4(i,2) = R*sin(theta) + L1*sinGamma;
    
    
    %Sistema Linear para descobrir W3 e W4
    A = [1 0 L1*sinGamma 0; 0 1 L1*cosGamma 0; 1 0 0 L3*sinAlpha; 0 1 0 L3*cosAlpha];
   
    B = [VP3(i,1); VP3(i,2); 0; 0];
    
    jk = A\B;
    
    VP4(i,1) = jk(1);
    VP4(i,2) = jk(2);
    
    W3(i,1) = jk(3);
    W4(i,1) = jk(4);
    
    %Cálculo da velocidade no ponto P6
    VP6(i,1) = W4(i,1)*L4*sinAlpha;
    VP6(i,2) = W4(i,1)*L4*cosAlpha;
    
    sinGamma2 = sinGamma+sin(15);
    cosGamma2 = cosGamma+sin(15);
    
    %Cálculo da velocidade no ponto P5
    VP5(i,1) = VP4(i,1)-W3(i,1)*((L1*sinGamma-L2)+(-L2*sinGamma2));
    VP5(i,2) = VP4(i,2)-W3(i,1)*((L1*sinGamma)+(-L2*cosGamma2));
    

    %Sistema Linear para descobrir velocidade no ponto P7, W5 e W6
    C = [1 0 L4*sinAlpha 0; 0 1 L4*cosAlpha 0; 1 0 0 L6*sinGamma; 0 1 0 L6*cosGamma];
    
    D = [VP5(i,1); VP5(i,2); VP6(i,1); VP6(i,2)];
    
    jk_2 = C\D;
    
    VP7(i,1) = jk_2(1);
    VP7(i,2) = jk_2(2);
    
    W5(i,1) = jk_2(3);
    W6(i,1) = jk_2(4);
    
    
    %Cálculo da velocidade no ponto P9
    VP9(i, 1) = VP5(i,1)+W5(i,1)*(L5*sinAlpha);
    VP9(i, 2) = VP5(i,2)+W5(i,1)*(L5*cosAlpha);
    
    
    Psi = 40;
    
    %Sistema Linear para decobrir velocidade no ponto P8, W7 e W8
    E = [1 0 -L6*sinAlpha 0; 0 1 -L6*cosAlpha 0; 1 0 0 d*sin(Psi); 0 1 0 d*cos(Psi)];
    
    F = [VP6(i,1); VP6(i,2); VP9(i,1); VP9(i,2)];
    
    jk_3 = E\F;
    
    VP8(i,1) = jk_3(1);
    VP8(i,2) = jk_3(2);
    
    W7(i,1) = jk_3(3);
    W8(i,1) = jk_3(4);

   
    

end

%Cálculo das posições do ponto P6
[T, PosP6] = ode45(@(t,PosP6)odeFunVP6(t, PosP6, w2,r2, L1, L2, R, x, L3, L4),tspan,rB0);

%Cálculo das posições do ponto P5
[T, PosP5] = ode45(@(t,PosP5)odeFunVP5(t, PosP5, w2,r2, L1, L2, R, x, L3, L4),tspan,rB0);

%Cálculo das posições do ponto P7
[T, PosP7] = ode45(@(t,PosP7)odeFunVP7(t, PosP7, w2,r2, L1, L2, R, x, L3, L4),tspan,rB0);

%Cálculo das posições do ponto P8
[T, PosP8] = ode45(@(t,PosP8)odeFunVP8(t, PosP8, w2,r2, L1, L2, R, x, L3, L4),tspan,rB0);

%Cálculo das posições do ponto P9
[T, PosP9] = ode45(@(t,PosP9)odeFunVP9(t, PosP9, w2,r2, L1, L2, R, x, L3, L4),tspan,rB0);



%Plots ---------------------------------------

figure (1)
plot(PosP3(:,1),PosP3(:,2))
xlabel('Pos P3 x (m)')
ylabel('Pos P3 y (m)')
grid

figure (2)
subplot(1,4,1)
plot(PosP4(:,1),PosP4(:,2))
grid
xlabel('Pos P4 x (m)')
ylabel('Pos P4 y (m)')

subplot(1,4,2)
plot(VP4(:,1),VP4(:,2))
grid
xlabel('Vel P4 x (m)')
ylabel('Vel P4 y (m)')


subplot(1,4,3)
plot(Theta,dCA)
grid
xlabel('\theta ')
ylabel('d (m)')

subplot(1,4,4)
plot(Theta,dCB)
grid
xlabel('\theta ')
ylabel('d (m)')


figure (3)
subplot(1,2,1)
plot(VP6(:,1),VP6(:,2))
xlabel('Vel P6 x (m)')
ylabel('Vel P6 y (m)')
grid

subplot(1,2,2)
plot(Theta, VP6(:,1),'b',Theta, VP6(:,2),'r')
xlabel('Vel P6 Theta x (m)')
ylabel('Vel P6 Theta y (m)')
grid


figure (4)

plot(PosP6(:,1),PosP6(:,2))
xlabel('Pos P6 x (m)')
ylabel('Pos P6 y (m)')
axis([-4 4 -4 4])
grid




figure (5)
plot(VP5(:,1),VP5(:,2))
xlabel('Vel P5 x (m)')
ylabel('Vel P5 y (m)')
grid

figure (6)
plot(PosP5(:,1),PosP5(:,2))
xlabel('Pos P5 x (m)')
ylabel('Pos P5 y (m)')
grid


figure (7)
plot(VP7(:,1),VP7(:,2))
xlabel('Vel P7 x (m)')
ylabel('Vel P7 y (m)')
grid

figure (8)
plot(VP8(:,1),VP8(:,2))
xlabel('Vel P8 x (m)')
ylabel('Vel P8 y (m)')
grid

figure (9)
plot(PosP7(:,1),PosP7(:,2))
xlabel('Pos P7 x (m)')
ylabel('Pos P7 y (m)')
axis([-90 50 -90 50])
grid

figure (10)
plot(PosP8(:,1),PosP8(:,2))
xlabel('Pos P8 x (m)')
ylabel('Pos P8 y (m)')
axis([-20 140 -190 40])
grid

figure (11)
plot(PosP9(:,1),PosP9(:,2))
xlabel('Pos P9 x (m)')
ylabel('Pos P9 y (m)')
axis([-120 60 -120 60])
grid

figure (12)
plot(PosP4(:,1),PosP4(:,2))
xlabel('Pos P4 x (m)')
ylabel('Pos P4 y (m)')
grid


%--------------------------------------------------------------------------


%Função que armazena as componentes da velocidade no ponto P3
function vA = calculaVelocidade(theta, w2, r2)

senot = sin(theta);
cost = cos(theta);

vA = [ -w2*r2*senot w2*r2*cost ];

%Função que irá integrar a velocidade do ponto P3 utilizando o ode45
function vB = odeFun(t,PosP3,w2,r2)

theta = w2*t;
vB = calculaVelocidade(theta, w2, r2);

vB = vB';

%Função que irá calcular junto ao Fsolve o sistema não linear
function F = funcaoMecanismoDeLee(x,theta,L1,L2,R)

xa = R*cos(theta);
ya = R*sin(theta);
xb = 4;
yb = 4;

F(1) = L1*sqrt(1-x(1)^2) - L2*sqrt(1-x(2)^2) - xa + xb;
F(2) = L1*x(1) - L2*x(2) - yb + ya;


%Função que armazena as componentes da velocidade no ponto P6
function vP6 = calculaVelocidade_VP6(theta, w2,r2, L1, L2, R, x, L3, L4)
X0 = [0 0];

senot = sin(theta);
cost = cos(theta);

VP3_x = -w2*r2*senot;
VP3_y = w2*r2*cost;

[x, ~] = fsolve(@(x)funcaoMecanismoDeLee(x,theta,L1,L2,R),X0);

sinGamma = x(1);
sinAlpha = x(2);
cosGamma = sqrt(1-sinGamma^2);
cosAlpha = sqrt(1-sinAlpha^2);

A = [1 0 L1*sinGamma 0; 0 1 L1*cosGamma 0; 1 0 0 L3*sinAlpha; 0 1 0 L3*cosAlpha];
    
B = [VP3_x; VP3_y; 0; 0];
    
jk = A\B;

VP4_x = jk(1);
VP4_y = jk(2);

W3_sistema = jk(3);
W4_sistema = jk(4);

sinGamma2 = sinGamma+sin(15);
cosGamma2 = cosGamma+sin(15);

vP6 = [ W4_sistema*L4*cosAlpha W4_sistema*L4*sinAlpha ];


%Função que irá integrar a velocidade do ponto P6 utilizando o ode45
function vBP6 = odeFunVP6(t, PosP6, w2,r2, L1, L2, R, x, L3, L4)

theta = w2*t;
vBP6 = calculaVelocidade_VP6(theta, w2,r2, L1, L2, R, x, L3, L4);

vBP6 = vBP6';

%Função que armazena as componentes da velocidade no ponto P5
function vP5 = calculaVelocidade_VP5(theta, w2,r2, L1, L2, R, x, L3, L4)

X0 = [0 0];

senot = sin(theta);
cost = cos(theta);

VP3_x = -w2*r2*senot;
VP3_y = w2*r2*cost;

[x, ~] = fsolve(@(x)funcaoMecanismoDeLee(x,theta,L1,L2,R),X0);

sinGamma = x(1);
sinAlpha = x(2);
cosGamma = sqrt(1-sinGamma^2);
cosAlpha = sqrt(1-sinAlpha^2);

A = [1 0 L1*sinGamma 0; 0 1 L1*cosGamma 0; 1 0 0 L3*sinAlpha; 0 1 0 L3*cosAlpha];
    
B = [VP3_x; VP3_y; 0; 0];
    
jk = A\B;

VP4_x = jk(1);
VP4_y = jk(2);

W3_sistema = jk(3);
W4_sistema = jk(4);

sinGamma2 = sinGamma+sin(15);
cosGamma2 = cosGamma+sin(15);

vP5 = [-(VP4_x-W3_sistema*((L1*sinGamma-L2)+(-L2*sinGamma2))) VP4_y-W3_sistema*((L1*sinGamma)+(-L2*cosGamma2))];

%Função que irá integrar a velocidade do ponto P5 utilizando o ode45
function vBP5 = odeFunVP5(t, PosP5, w2,r2, L1, L2, R, x, L3, L4)

theta = w2*t;
vBP5 = calculaVelocidade_VP5(theta, w2,r2, L1, L2, R, x, L3, L4);

vBP5 = vBP5';




%Função que armazena as componentes da velocidade no ponto P7
function vP7 = calculaVelocidade_VP7(theta, w2,r2, L1, L2, R, x, L3, L4)

X0 = [0 0];

senot = sin(theta);
cost = cos(theta);

VP3_x = -w2*r2*senot;
VP3_y = w2*r2*cost;

[x, ~] = fsolve(@(x)funcaoMecanismoDeLee(x,theta,L1,L2,R),X0);

sinGamma = x(1);
sinAlpha = x(2);
cosGamma = sqrt(1-sinGamma^2);
cosAlpha = sqrt(1-sinAlpha^2);

A = [1 0 L1*sinGamma 0; 0 1 L1*cosGamma 0; 1 0 0 L3*sinAlpha; 0 1 0 L3*cosAlpha];
   
B = [VP3_x; VP3_y; 0; 0];
    
jk = A\B;

VP4_x = jk(1);
VP4_y = jk(2);

W3_sistema = jk(3);
W4_sistema = jk(4);

sinGamma2 = sinGamma+sin(15);
cosGamma2 = cosGamma+sin(15);

VP6_x = W4_sistema*L4*sinAlpha;
VP6_y = W4_sistema*L4*cosAlpha;

Psi = 40;
L5 = 29.2;
L6 = 2.7;
d = 2;


VP5_x = VP4_x-W3_sistema*((L1*sinGamma-L2)+(-L2*sinGamma2));
VP5_y = VP4_y-W3_sistema*((L1*sinGamma)+(-L2*cosGamma2));


C = [1 0 L4*sinAlpha 0; 0 1 L4*cosAlpha 0; 1 0 0 L6*sinGamma; 0 1 0 L6*cosGamma];

D = [VP5_x; VP5_y; VP6_x; VP6_y];

jk_2 = C\D;

VP7_x = jk_2(1);
VP7_y = jk_2(2);

W5_sistema = jk_2(3);
W6_sistema = jk_2(4);

VP9_x = VP5_x+W5_sistema*(L5*sinAlpha);
VP9_y = VP5_y+W5_sistema*(L5*cosAlpha);

    
E = [1 0 -L6*sinAlpha 0; 0 1 -L6*cosAlpha 0; 1 0 0 d*sin(Psi); 0 1 0 d*cos(Psi)];

F = [VP6_x; VP6_y; VP9_x; VP9_y];

jk_3 = E\F;

VP8_x = jk_3(1);
VP8_y = jk_3(2);

W7_sistema = jk_3(3);
W8_sistema = jk_3(4);


vP7 = [VP5_x-(L4*W5_sistema*sinAlpha) VP5_y-(L4*W5_sistema*cosAlpha)];


%Função que irá integrar a velocidade do ponto P7 utilizando o ode45
function vBP7 = odeFunVP7(t, PosP7, w2,r2, L1, L2, R, x, L3, L4)

theta = w2*t;
vBP7 = calculaVelocidade_VP7(theta, w2,r2, L1, L2, R, x, L3, L4);

vBP7 = vBP7';






%Função que armazena as componentes da velocidade no ponto P8
function vP8 = calculaVelocidade_VP8(theta, w2,r2, L1, L2, R, x, L3, L4)

X0 = [0 0];

senot = sin(theta);
cost = cos(theta);

VP3_x = -w2*r2*senot;
VP3_y = w2*r2*cost;

[x, ~] = fsolve(@(x)funcaoMecanismoDeLee(x,theta,L1,L2,R),X0);

sinGamma = x(1);
sinAlpha = x(2);
cosGamma = sqrt(1-sinGamma^2);
cosAlpha = sqrt(1-sinAlpha^2);

A = [1 0 L1*sinGamma 0; 0 1 L1*cosGamma 0; 1 0 0 L3*sinAlpha; 0 1 0 L3*cosAlpha];
   
B = [VP3_x; VP3_y; 0; 0];
    
jk = A\B;

VP4_x = jk(1);
VP4_y = jk(2);

W3_sistema = jk(3);
W4_sistema = jk(4);

sinGamma2 = sinGamma+sin(15);
cosGamma2 = cosGamma+sin(15);

VP6_x = W4_sistema*L4*sinAlpha;
VP6_y = W4_sistema*L4*cosAlpha;

Psi = 40;
L5 = 29.2;
L6 = 2.7;
d = 2;


VP5_x = VP4_x-W3_sistema*((L1*sinGamma-L2)+(-L2*sinGamma2));
VP5_y = VP4_y-W3_sistema*((L1*sinGamma)+(-L2*cosGamma2));

C = [1 0 -L4*sinAlpha 0; 0 1 -L4*cosAlpha 0; 1 0 0 L6*sinGamma; 0 1 0 L6*cosGamma];

D = [VP5_x; VP5_y; VP6_x; VP6_y];

jk_2 = C\D;

VP7_x = jk_2(1);
VP7_y = jk_2(2);

W5_sistema = jk_2(3);
W6_sistema = jk_2(4);

VP9_x = VP5_x+W5_sistema*(L5*sinAlpha);
VP9_y = VP5_y+W5_sistema*(L5*cosAlpha);

    
E = [1 0 -L6*sinAlpha 0; 0 1 -L6*cosAlpha 0; 1 0 0 d*sin(Psi); 0 1 0 d*cos(Psi)];

F = [VP6_x; VP6_y; VP9_x; VP9_y];

jk_3 = E\F;

VP8_x = jk_3(1);
VP8_y = jk_3(2);

W7_sistema = jk_3(3);
W8_sistema = jk_3(4);

vP8 = [VP6_x+W7_sistema*L6*sinAlpha VP6_y+W7_sistema*L6*cosAlpha];


%Função que irá integrar a velocidade do ponto P8 utilizando o ode45
function vBP8 = odeFunVP8(t, PosP8, w2,r2, L1, L2, R, x, L3, L4)

theta = w2*t;
vBP8 = calculaVelocidade_VP8(theta, w2,r2, L1, L2, R, x, L3, L4);

vBP8 = vBP8';


%Função que armazena as componentes da velocidade no ponto P9
function vP9 = calculaVelocidade_VP9(theta, w2,r2, L1, L2, R, x, L3, L4)

X0 = [0 0];

senot = sin(theta);
cost = cos(theta);

VP3_x = -w2*r2*senot;
VP3_y = w2*r2*cost;

[x, ~] = fsolve(@(x)funcaoMecanismoDeLee(x,theta,L1,L2,R),X0);

sinGamma = x(1);
sinAlpha = x(2);
cosGamma = sqrt(1-sinGamma^2);
cosAlpha = sqrt(1-sinAlpha^2);

A = [1 0 L1*sinGamma 0; 0 1 L1*cosGamma 0; 1 0 0 L3*sinAlpha; 0 1 0 L3*cosAlpha];
   
B = [VP3_x; VP3_y; 0; 0];
    
jk = A\B;

VP4_x = jk(1);
VP4_y = jk(2);

W3_sistema = jk(3);
W4_sistema = jk(4);

sinGamma2 = sinGamma+sin(15);
cosGamma2 = cosGamma+sin(15);

VP6_x = W4_sistema*L4*sinAlpha;
VP6_y = W4_sistema*L4*cosAlpha;

Psi = 40;
L5 = 29.2;
L6 = 2.7;
d = 2;

VP5_x = VP4_x-W3_sistema*((L1*sinGamma-L2)+(-L2*sinGamma2));
VP5_y = VP4_y-W3_sistema*((L1*sinGamma)+(-L2*cosGamma2));

C = [1 0 -L4*sinAlpha 0; 0 1 -L4*cosAlpha 0; 1 0 0 L6*sinGamma; 0 1 0 L6*cosGamma];

D = [VP5_x; VP5_y; VP6_x; VP6_y];

jk_2 = C\D;

VP7_x = jk_2(1);
VP7_y = jk_2(2);

W5_sistema = jk_2(3);
W6_sistema = jk_2(4);

VP9_x = VP5_x+W5_sistema*(L5*sinAlpha);
VP9_y = VP5_y+W5_sistema*(L5*cosAlpha);


E = [1 0 -L6*sinAlpha 0; 0 1 -L6*cosAlpha 0; 1 0 0 d*sin(Psi); 0 1 0 d*cos(Psi)];

F = [VP6_x; VP6_y; VP9_x; VP9_y];

jk_3 = E\F;

VP8_x = jk_3(1);
VP8_y = jk_3(2);

W7_sistema = jk_3(3);
W8_sistema = jk_3(4);

vP9 = [VP8_x+W8_sistema*d*sin(Psi) VP8_y+W8_sistema*d*cos(Psi)];







%Função que irá integrar a velocidade do ponto P9 utilizando o ode45
function vBP9 = odeFunVP9(t, PosP9, w2,r2, L1, L2, R, x, L3, L4)

theta = w2*t;
vBP9 = calculaVelocidade_VP9(theta, w2,r2, L1, L2, R, x, L3, L4);

vBP9 = vBP9';


