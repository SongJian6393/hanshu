
function [f,df] = hanshu2_0(XDENS)

% XDENS =  [50000,10000,1000,100];

DataC_0;

Xigma_YS0 = 180;

Xigma_sat = XDENS(1);

Beta = XDENS(2);

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

f = sum([Xigma_2 - Data(51:142,2); Xigma_3 - Data(150:241,2)].^2);

dlta = 0.1;

for i = 1:2
    
    XXDENS = XDENS;
    
    XXDENS(i) = XXDENS(i)+dlta;
    
    Xigma_sat = XXDENS(1);

    Beta = XXDENS(2);
    

    alpha_11 = C1/gama1*(1-exp(-gama1*evpe_1));

    alpha_12 = C2/gama2*(1-exp(-gama2*evpe_1));
    

    alpha_21 = (alpha_11(end)+C1/gama1)*(exp(-gama1*(evp_1(end)-evp_2))) - C1/gama1;

    alpha_22 = (alpha_12(end)+C2/gama2)*(exp(-gama2*(evp_1(end)-evp_2))) - C2/gama2;

    alpha_2 = alpha_21 + alpha_22;
    
    Xigma_2 = alpha_2 - Xigma_YS_2;


    alpha_31 = (alpha_21(end)-C1/gama1)*(exp(-gama1*(evp_3-evp_2(end)))) + C1/gama1;

    alpha_32 = (alpha_22(end)-C2/gama2)*(exp(-gama2*(evp_3-evp_2(end)))) + C2/gama2;

    alpha_3 = alpha_31 + alpha_32;

    Xigma_3 = alpha_3 + Xigma_YS_3;   
   

    f_d = sum([Xigma_2 - Data(51:142,2); Xigma_3 - Data(150:241,2)].^2)

    df(i) = (f_d -f)/dlta;
   
end
 
