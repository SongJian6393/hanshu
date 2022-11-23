clear all

OBJE=@hanshu1_0;

N = 4;

LB =  [10000,5000,100,1];
UB =  [150000,50000,5000,1000];

x0 =  [50000,10000,1000,100];

OP=optimset( 'Algorithm','trust-region-reflective','GradObj','on','disp','iter','MaxFunEvals',200,'TolFun',1e-20,'PlotFcns',@optimplotfval);

[XOPT] = fmincon(OBJE,x0,[],[],[],[],LB,UB,[],OP)

figure,

DataC_0;

hold on

% plot(Data(1:50,1),Data(1:50,2))
% plot(Data(50:140,1),Data(50:140,2))
% plot(Data(149:239,1),Data(149:239,2))

Xigma_YS0 = 180;

Xigma_sat = 300;

Beta = 2;

C1 = XOPT(1); C2 = XOPT(2); 

gama1 = XOPT(3); gama2 = XOPT(4);
 
evp_1 = Data(1:43,1);

evpe_1 = evp_1 ;

Xigma_YS_1 = Xigma_YS0 + Xigma_sat*(1-exp(-Beta*evpe_1));        %voce 模型

alpha_11 = C1/gama1*(1-exp(-gama1*evpe_1));

alpha_12 = C2/gama2*(1-exp(-gama2*evpe_1));

alpha_1 = alpha_11 + alpha_12;    % chaboche模型

Xigma_1 = alpha_1 + Xigma_YS_1;         % 总应力

evp_2 = Data(51:142,1);

evpe_2 = evpe_1(end) + sqrt((evp_2 - evp_2(1)).^2);

Xigma_YS_2 = Xigma_YS0 + Xigma_sat*(1-exp(-Beta*evpe_2));

alpha_21 = (alpha_11(end)+C1/gama1)*(exp(-gama1*(evp_1(end)-evp_2))) - C1/gama1;

alpha_22 = (alpha_12(end)+C2/gama2)*(exp(-gama2*(evp_1(end)-evp_2))) - C2/gama2;

alpha_2 = alpha_21 + alpha_22;

Xigma_2 = alpha_2 - Xigma_YS_2;

evp_3 = Data(150:241,1);

evpe_3 = evpe_2(end) + (evp_3 - evp_3(1));

Xigma_YS_3 = Xigma_YS0 + Xigma_sat*(1-exp(-Beta*evpe_3));

alpha_31 = (alpha_21(end)-C1/gama1)*(exp(-gama1*(evp_3-evp_2(end)))) + C1/gama1;

alpha_32 = (alpha_22(end)-C2/gama2)*(exp(-gama2*(evp_3-evp_2(end)))) + C2/gama2;

alpha_3 = alpha_31 + alpha_32;

Xigma_3 = alpha_3 + Xigma_YS_3; 

f = sum([Xigma_2 - Data(51:142,2); Xigma_3 - Data(150:241,2)].^2)


hold on 
plot (evp_1,Xigma_1,'r--');
plot (evp_2,Xigma_2,'r--');
plot (evp_3,Xigma_3,'r--');






