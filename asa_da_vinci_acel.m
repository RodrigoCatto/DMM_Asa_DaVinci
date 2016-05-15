function asa_da_vinci_acel(w1)

r1 = 1.35;
r2 = 3.9;

w2 = (-w1*r1)/r2;
v2 = r1*w2;


R = sqrt(2);
L1 = sqrt(20);
L2 = sqrt(66);
L3 = 8.3;
L4 = 18.1;
L5 = 29.2;
L6 = 2.7;
L7 = L6;
d = 2;


numPontos = 100;
Theta = linspace(0,360,numPontos);
X1 = zeros(numPontos,1);
X2 = zeros(numPontos,1);
VP3 = zeros(numPontos,2);
VelP4 = zeros(numPontos,2);
VP6 = zeros(numPontos,2);
VP5 = zeros(numPontos,2);
VP7 = zeros(numPontos,2);
VP8 = zeros(numPontos,2);
VP9 = zeros(numPontos,2);

AP3 = zeros(numPontos,2);
AP4 = zeros(numPontos,2);
AP5 = zeros(numPontos,2);
AP6 = zeros(numPontos,2);
AP7 = zeros(numPontos,2);
AP8 = zeros(numPontos,2);
AP9 = zeros(numPontos,2);
VP3 = zeros(numPontos,2);

Forca = zeros(numPontos,2);
Far = zeros(numPontos,1);
M2 = zeros(numPontos,1);
M3 = zeros(numPontos,1);


W3 = zeros(numPontos,1);
W4 = zeros(numPontos,1);
W5 = zeros(numPontos,1);
W6 = zeros(numPontos,1);
W7 = zeros(numPontos,1);
W8 = zeros(numPontos,1);


W3_ponto = zeros(numPontos,1);
W4_ponto = zeros(numPontos,1);
W5_ponto = zeros(numPontos,1);
W6_ponto = zeros(numPontos,1);
W7_ponto = zeros(numPontos,1);
W8_ponto = zeros(numPontos,1);

dCA = zeros(numPontos,1);
dCB = zeros(numPontos,1);

f = abs(w2/2/pi);
tspan = [ 0 1/f ];
rB0 = [ r2 0 ];

[T, PosP3] = ode45(@(t,PosP3)odeFun(t,PosP3,w2,r2),tspan,rB0);


for i = 1: numPontos
    theta = Theta(i)*pi/180;
    
    senot = sin(theta);
    cost = cos(theta);

    VP3(i,1) = -w2*r2*senot
    VP3(i,2) = w2*r2*cost
    
    if i == 1
        X0 = [0 0];
    else
        X0 = [X1(i-1) X2(i-1)];
    end
    [x, ~] = fsolve(@(x)funcaoMecanismoDeLee(x,theta,L1,L2,R),X0);
    X1(i) = x(1);
    X2(i) = x(2);
    sinGamma = x(1);
    sinAlpha = x(2);
    cosGamma = sqrt(1-sinGamma^2);
    cosAlpha = sqrt(1-sinAlpha^2);
    PosP4(i,1) = R*cos(theta) - L1*cosGamma;
    PosP4(i,2) = R*sin(theta) + L1*sinGamma;
    
    
    
    A = [1 0 L1*sinGamma 0; 0 1 L1*cosGamma 0; 1 0 0 L3*sinAlpha; 0 1 0 L3*cosAlpha];
    %x = [VP4x; VP4y; W3; W4];
    B = [VP3(i,1); VP3(i,2); 0; 0];
    
    jk = A\B;
    
    VP4(i,1) = jk(1);
    VP4(i,2) = jk(2);
    
    W3(i,1) = jk(3);
    W4(i,1) = jk(4);
    
    VP6(i,1) = W4(i,1)*L4*sinAlpha;
    VP6(i,2) = W4(i,1)*L4*cosAlpha;
    
    sinGamma2 = sinGamma+sin(15);
    cosGamma2 = cosGamma+sin(15);
    
    VP5(i,1) = VP4(i,1)-W3(i,1)*((L1*sinGamma-L2)+(-L2*sinGamma2));
    VP5(i,2) = VP4(i,2)-W3(i,1)*((L1*sinGamma)+(-L2*cosGamma2));
    
    dCA(i) = sqrt((PosP4(i,1)-R*cos(theta))^2 + (PosP4(i,2)-R*sin(theta))^2);
    xb = 4;
    yb = 4;
   
    dCB(i) = sqrt((PosP4(i,1)-xb)^2 + (PosP4(i,2)-yb)^2);

    
    C = [1 0 L4*sinAlpha 0; 0 1 L4*cosAlpha 0; 1 0 0 L6*sinGamma; 0 1 0 L6*cosGamma];
    %x = [VP4x; VP4y; W3; W4];
    D = [VP5(i,1); VP5(i,2); VP6(i,1); VP6(i,2)];
    
    jk_2 = C\D;
    
    VP7(i,1) = jk_2(1);
    VP7(i,2) = jk_2(2);
    
    W5(i,1) = jk_2(3);
    W6(i,1) = jk_2(4);
    
    VP9(i, 1) = VP5(i,1)+W5(i,1)*(L5*sinAlpha);
    VP9(i, 2) = VP5(i,2)+W5(i,1)*(L5*cosAlpha);
    
    
    Psi = 40;
    
    E = [1 0 -L6*sinAlpha 0; 0 1 -L6*cosAlpha 0; 1 0 0 d*sin(Psi); 0 1 0 d*cos(Psi)];
    %x = [VP4x; VP4y; W3; W4];
    F = [VP6(i,1); VP6(i,2); VP9(i,1); VP9(i,2)];
    
    jk_3 = E\F;
    
    VP8(i,1) = jk_3(1);
    VP8(i,2) = jk_3(2);
    
    W7(i,1) = jk_3(3);
    W8(i,1) = jk_3(4);
    
    vetor_dist_pontos_P3 = [r2*cost r2*senot 0];
    vel_angular_P3 = [0 0 w2];

    acel_pt_2 = cross(vel_angular_P3, cross(vetor_dist_pontos_P3, vel_angular_P3));


    AP3(i,1) = [acel_pt_2(1,1)]; 
    AP3(i,2) = [acel_pt_2(1,2)];
    
    G = [1 0 L1*sinGamma 0;1 0 0 L3*sinAlpha;0 1 L1*cosGamma 0;0 1 0 L3*cosAlpha]
    H = [AP3(i,1)+(W3(i,1)^2)*L1*cosGamma; (W4(i,1)^2)*L3*cosAlpha; AP3(i,2)-(W3(i,1)^2)*L1*sinGamma; (W3(i,1)^2)*L3*cosAlpha];
    
    jk_4 = G\H;
    
    

    AP4(i,1) = jk_4(1);
    AP4(i,2) = jk_4(2);
    
    W3_ponto(i,1) = jk_4(3);
    W4_ponto(i,1) = jk_4(4);
    
    rp6_rp0 = [-L4*cosAlpha L4*sinAlpha 0];
    Cross_P6_P0 = cross([0 0 W4_ponto(i,1)],rp6_rp0)+ cross([0 0 W4(i,1)], cross([0 0 W4(i,1)], rp6_rp0)); 
    
    AP6(i,1) = Cross_P6_P0(1,1);
    AP6(i,2) = Cross_P6_P0(1,2);
    
    rp5_rp4 = [((L1*sinGamma-L2)+(-L2*sinGamma2)) ((L1*sinGamma)+(-L2*cosGamma2)) 0];
    Cross_P5_P4 = cross([0 0 W3_ponto(i,1)], rp5_rp4) + cross([0 0 W3(i,1)], cross([0 0 W3(i,1)], rp5_rp4));
    
    AP5(i,1) = AP4(i,1)+ Cross_P5_P4(1,1);
    AP5(i,2) = AP4(i,2)+ Cross_P5_P4(1,2);
    
    I = [1 0 L4*sinAlpha 0;1 0 0 L6*sinGamma;0 1 L4*cosAlpha 0;0 1 0 L6*cosGamma]
    J = [AP5(i,1)+((W5(i,1))^2)*L4*cosAlpha; AP6(i,1)+(W6(i,1)^2)*L6*cosGamma; AP5(i,2)-(W5(i,1)^2)*L4*sinAlpha; AP6(i,2)- (W6(i,1)^2)*L6*sinGamma];
    
    jk_5 = I\J;
    
    AP7(i,1) = jk_5(1);
    AP7(i,2) = jk_5(2);
    
    W5_ponto(i,1) = jk_5(3);
    W6_ponto(i,1) = jk_5(4);
    
    rp9_rp5 = [L5*cosAlpha L5*sinAlpha 0]
    Cross_P9_P5 = cross([0 0 W5_ponto(i,1)], rp9_rp5) + cross([0 0 W5(i,1)], cross([0 0 W5(i,1)], rp9_rp5));
    
    AP9(i,1) = AP5(i,1) + Cross_P9_P5(1,1);
    AP9(i,2) = AP5(i,2) + Cross_P9_P5(1,2);
    
    G = [1 0 -L7*sinAlpha 0;1 0 0 -d*sin(Psi);0 1 -L7*cosAlpha 0;0 1 0 d*cos(Psi)]
    H = [AP6(i,1)- (W7(i,1)^2)*L7*cosAlpha; AP9(i,1)+(W8(i,1)^2)*d*cos(Psi); AP6(i,2)+(W7(i,1)^2)*L7*sinAlpha; AP9(i,2)+(W8(i,1)^2)*d*sin(Psi)];
    
    jk_6 = G\H;
    
    AP8(i,1) = jk_6(1);
    AP8(i,2) = jk_6(2);
    
    W7_ponto(i,1) = jk_6(3);
    W8_ponto(i,1) = jk_6(4);
    
    Peso = 100;
    Massa_barra = 10;
    Pot = 90;
    
    
    Far(i,1) = (((Massa_barra*(L4+L5)^2)/12)*W4_ponto(i,1) + (L4+L5)*Peso*cosAlpha - L3*(Pot/v2))/(L5+L4);
    
    
end



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




function vA = calculaVelocidade(theta, w2, r2)

senot = sin(theta);
cost = cos(theta);

vA = [ -w2*r2*senot w2*r2*cost ];


function vB = odeFun(t,PosP3,w2,r2)

theta = w2*t;
vB = calculaVelocidade(theta, w2, r2);

vB = vB';


function F = funcaoMecanismoDeLee(x,theta,L1,L2,R)

xa = R*cos(theta);
ya = R*sin(theta);
xb = 4;
yb = 4;

F(1) = L1*sqrt(1-x(1)^2) - L2*sqrt(1-x(2)^2) - xa + xb;
F(2) = L1*x(1) - L2*x(2) - yb + ya;


