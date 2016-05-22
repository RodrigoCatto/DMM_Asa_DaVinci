function asa_da_vinci_acel(w1)

r1 = 1.35;
r2 = 3.9;

w2 = (-w1*r1)/r2; %Vel angular engrenagem maior

v2 = r1*w2;

%Tamanho das barras do mecanismo
R = sqrt(2);
L1 = sqrt(20);
L2 = sqrt(66);
L3 = 8.3;
L4 = 18.1;
L5 = 29.2;
L6 = 2.7;
L7 = L6;
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

%Matrizes que armazenam os valores das acelerações
AP3 = zeros(numPontos,2);
AP4 = zeros(numPontos,2);
AP5 = zeros(numPontos,2);
AP6 = zeros(numPontos,2);
AP7 = zeros(numPontos,2);
AP8 = zeros(numPontos,2);
AP9 = zeros(numPontos,2);
VP3 = zeros(numPontos,2);

%Matrizes que armazenam os valores das forças
Forca = zeros(numPontos,2);
Far = zeros(numPontos,1);
M2 = zeros(numPontos,1);
M3 = zeros(numPontos,1);

%Matrizes que armazenam os valores das velocidades angulares
W3 = zeros(numPontos,1);
W4 = zeros(numPontos,1);
W5 = zeros(numPontos,1);
W6 = zeros(numPontos,1);
W7 = zeros(numPontos,1);
W8 = zeros(numPontos,1);

%Matrizes que armazenam os valores das acelerações angulares
W3_ponto = zeros(numPontos,1);
W4_ponto = zeros(numPontos,1);
W5_ponto = zeros(numPontos,1);
W6_ponto = zeros(numPontos,1);
W7_ponto = zeros(numPontos,1);
W8_ponto = zeros(numPontos,1);

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
    VP3(i,1) = -w2*r2*senot
    VP3(i,2) = w2*r2*cost
    
    %Cálculo do sistema não linear
    if i == 1
        X0 = [0 0];
    else
        X0 = [X1(i-1) X2(i-1)];
    end
    [x, ~] = fsolve(@(x)funcaoMecanismoDeLee(x,theta,L1,L2,R),X0);
    X1(i) = x(1);
    X2(i) = x(2);
    %Resultados do sistema não linear (sinGamma e sinAlpha)
    sinGamma = x(1);
    sinAlpha = x(2);
    cosGamma = sqrt(1-sinGamma^2);
    cosAlpha = sqrt(1-sinAlpha^2);
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
    
    
    vetor_dist_pontos_P3 = [r2*cost r2*senot 0];
    vel_angular_P3 = [0 0 w2];
    %Cálculo do cross entre os vetores para a aceleração do ponto P3
    acel_pt_2 = cross(vel_angular_P3, cross(vetor_dist_pontos_P3, vel_angular_P3));

    %Cálculo da aceleração do ponto P3
    AP3(i,1) = [acel_pt_2(1,1)]; 
    AP3(i,2) = [acel_pt_2(1,2)];
    
    %Sistema Linear para decobrir aceleração no ponto P4, aceleração
    %angular do corpo 3 e do corpo 4
    G = [1 0 L1*sinGamma 0;1 0 0 L3*sinAlpha;0 1 L1*cosGamma 0;0 1 0 L3*cosAlpha]
    H = [AP3(i,1)+(W3(i,1)^2)*L1*cosGamma; (W4(i,1)^2)*L3*cosAlpha; AP3(i,2)-(W3(i,1)^2)*L1*sinGamma; (W3(i,1)^2)*L3*cosAlpha];
    
    jk_4 = G\H;

    AP4(i,1) = jk_4(1);
    AP4(i,2) = jk_4(2);
    
    W3_ponto(i,1) = jk_4(3);
    W4_ponto(i,1) = jk_4(4);
    
    rp6_rp0 = [-L4*cosAlpha L4*sinAlpha 0];
    %Cálculo do cross entre os vetores para a aceleração do ponto P6
    Cross_P6_P0 = cross([0 0 W4_ponto(i,1)],rp6_rp0)+ cross([0 0 W4(i,1)], cross([0 0 W4(i,1)], rp6_rp0)); 
    
    %Cálculo da aceleração do ponto P6
    AP6(i,1) = Cross_P6_P0(1,1);
    AP6(i,2) = Cross_P6_P0(1,2);
    
    
    
    rp5_rp4 = [((L1*sinGamma-L2)+(-L2*sinGamma2)) ((L1*sinGamma)+(-L2*cosGamma2)) 0];
    %Cálculo do cross entre os vetores para a aceleração do ponto P5
    Cross_P5_P4 = cross([0 0 W3_ponto(i,1)], rp5_rp4) + cross([0 0 W3(i,1)], cross([0 0 W3(i,1)], rp5_rp4));
    
    %Cálculo da aceleração do ponto P5
    AP5(i,1) = AP4(i,1)+ Cross_P5_P4(1,1);
    AP5(i,2) = AP4(i,2)+ Cross_P5_P4(1,2);
    
    
    %Sistema Linear para decobrir aceleração no ponto P7, aceleração
    %angular do corpo 5 e do corpo 6
    I = [1 0 L4*sinAlpha 0;1 0 0 L6*sinGamma;0 1 L4*cosAlpha 0;0 1 0 L6*cosGamma]
    J = [AP5(i,1)+((W5(i,1))^2)*L4*cosAlpha; AP6(i,1)+(W6(i,1)^2)*L6*cosGamma; AP5(i,2)-(W5(i,1)^2)*L4*sinAlpha; AP6(i,2)- (W6(i,1)^2)*L6*sinGamma];
    
    jk_5 = I\J;
    
    AP7(i,1) = jk_5(1);
    AP7(i,2) = jk_5(2);
    
    W5_ponto(i,1) = jk_5(3);
    W6_ponto(i,1) = jk_5(4);
    
    
    
    rp9_rp5 = [L5*cosAlpha L5*sinAlpha 0]
    %Cálculo do cross entre os vetores para a aceleração do ponto P9
    Cross_P9_P5 = cross([0 0 W5_ponto(i,1)], rp9_rp5) + cross([0 0 W5(i,1)], cross([0 0 W5(i,1)], rp9_rp5));
    
    %Cálculo da aceleração do ponto P9
    AP9(i,1) = AP5(i,1) + Cross_P9_P5(1,1);
    AP9(i,2) = AP5(i,2) + Cross_P9_P5(1,2);
    
    
    %Sistema Linear para decobrir aceleração no ponto P8, aceleração
    %angular do corpo 7 e do corpo 8
    G = [1 0 -L7*sinAlpha 0;1 0 0 -d*sin(Psi);0 1 -L7*cosAlpha 0;0 1 0 d*cos(Psi)]
    H = [AP6(i,1)- (W7(i,1)^2)*L7*cosAlpha; AP9(i,1)+(W8(i,1)^2)*d*cos(Psi); AP6(i,2)+(W7(i,1)^2)*L7*sinAlpha; AP9(i,2)+(W8(i,1)^2)*d*sin(Psi)];
    
    jk_6 = G\H;
    
    AP8(i,1) = jk_6(1);
    AP8(i,2) = jk_6(2);
    
    W7_ponto(i,1) = jk_6(3);
    W8_ponto(i,1) = jk_6(4);
    
    %Cálculos para a dinâmica --------------------------------------------
    
    Peso = 100;
    Massa_barra = 10;
    Pot = 90;
   
    %Cálculo da Força do ar
    Far(i,1) = (((Massa_barra*(L4+L5)^2)/12)*W4_ponto(i,1) + (L4+L5)*Peso*cosAlpha - L3*(Pot/v2))/(L5+L4);
    
    
end


%Plots --------------------------------------------------------------------


figure (1)
plot(AP3(:,1),AP3(:,2))
title('Aceleração Ponto 3')
xlabel('Acel P3 x (m/s^2)')
ylabel('Acel P3 y (m/s^2)')
grid

figure (2)
plot(AP4(:,1),AP4(:,2))
title('Aceleração Ponto 4')
xlabel('Acel P4 x (m/s^2)')
ylabel('Acel P4 y (m/s^2)')
grid

figure (3)
plot(AP6(:,1),AP6(:,2))
title('Aceleração Ponto 6')
xlabel('Acel P6 x (m/s^2)')
ylabel('Acel P6 y (m/s^2)')
grid

figure (3)
plot(AP5(:,1),AP5(:,2))
title('Aceleração Ponto 5')
xlabel('Acel P5 x (m/s^2)')
ylabel('Acel P5 y (m/s^2)')
grid

figure (4)
plot(AP7(:,1),AP7(:,2))
title('Aceleração Ponto 7')
xlabel('Acel P7 x (m/s^2)')
ylabel('Acel P7 y (m/s^2)')
grid



figure (5)
plot(AP9(:,1),AP9(:,2))
title('Aceleração Ponto 9')
xlabel('Acel P9 x (m/s^2)')
ylabel('Acel P9 y (m/s^2)')
grid


figure (6)
plot(AP8(:,1),AP8(:,2))
title('Aceleração Ponto 8')
xlabel('Acel P8 x (m/s^2)')
ylabel('Acel P8 y (m/s^2)')
grid

figure (7)
plot(AP6(:,1),AP6(:,2))
title('Aceleração Ponto 6')
xlabel('Acel P6 x (m/s^2)')
ylabel('Acel P6 y (m/s^2)')
grid


figure (8)
plot(Theta,  Far(:,1))
title('Força Ar')
xlabel('Theta (Nm)')
ylabel('Força (N)')
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



