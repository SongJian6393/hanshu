clear all

OBJE=@hanshu2_0;

N = 2;

LB =  [10,1];
UB =  [1000,100];

x0 =  [319,1.95];

OP=optimset( 'Algorithm','trust-region-reflective','GradObj','on','disp','iter','MaxFunEvals',200,'TolFun',1e-20,'PlotFcns',@optimplotfval);

[XOPT] = fmincon(OBJE,x0,[],[],[],[],LB,UB,[],OP)

figure,

DataC_0;

hold on

% plot(Data(1:50,1),Data(1:50,2))
% plot(Data(53:140,1),Data(53:140,2))
% plot(Data(152:239,1),Data(152:239,2))

Xigma_YS0 = 180;

Xigma_sat = XOPT(1);

Beta = XOPT(2);

C1 = 110470; C2 = 26990; 

gama1 = 1130; gama2 = 220;

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


% hold on 
% plot (evp_1,Xigma_1,'*');
% plot (evp_2,Xigma_2,'*');
% plot (evp_3,Xigma_3,'*');
hold on 
plot (evp_1,Xigma_1,'r--');
plot (evp_2,Xigma_2,'r--');
plot (evp_3,Xigma_3,'r--');
legend('实验曲线','拟合曲线')
x_val=get(gca,'XTick');
x_str=num2str(x_val');
set(gca,'XTicklabel',x_str);
xlabel('塑性应变');
ylabel('应力/MPa');




