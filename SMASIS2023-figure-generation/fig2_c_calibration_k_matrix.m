clc;close all;clear;
cd ./data4plot/Kmatrix/
load('C_1.mat');
load('C_2.mat');
load('C_3.mat');
load('C_4.mat');
load('ness_temp1.mat');
load('ness_temp2.mat');
load('ness_temp3.mat');
load('ness_temp4.mat');
load('pwm1.mat');
load('pwm2.mat');
load('pwm3.mat');
load('pwm4.mat');
cd ../../

%%
indx1 = find(ness_temp1>50);
C_1 = C_1(indx1);
ness_temp1 = ness_temp1(indx1);
pwm1 = pwm1(indx1);
indx2 = find(ness_temp2>50);
C_2 = C_2(indx2);
ness_temp2 = ness_temp2(indx2);
pwm2 = pwm2(indx2);
indx3 = find(ness_temp3>50);
C_3 = C_3(indx3);
ness_temp3 = ness_temp3(indx3);
pwm3 = pwm3(indx3);
indx4 = find(ness_temp4>50);
C_4 = C_4(indx4);
ness_temp4 = ness_temp4(indx4);
pwm4 = pwm4(indx4);

%
cool_temp1 = ness_temp1(pwm1==0);
cool_C_1 = C_1(pwm1==0);
cool_mdl_1 = fitlm(cool_temp1,cool_C_1);

cool_mdl_1 = cool_mdl_1.Coefficients;
cool_int_1 = table2array(cool_mdl_1(1,1));
cool_slope_1 = table2array(cool_mdl_1(2,1));

lc1_x = 50:1:108;
lc1_y = cool_int_1 + cool_slope_1*lc1_x;
%
heat_temp1 = ness_temp1(pwm1~=0);
heat_C_1 = C_1(pwm1~=0);
heat_mdl_1 = fitlm(heat_temp1,heat_C_1);

heat_mdl_1 = heat_mdl_1.Coefficients;
heat_int_1 = table2array(heat_mdl_1(1,1));
heat_slope_1 = table2array(heat_mdl_1(2,1));

lh1_x = 50:1:108;
lh1_y = heat_int_1 + heat_slope_1*lh1_x;
%
cool_temp2 = ness_temp2(pwm2==0);
cool_C_2 = C_2(pwm2==0);
cool_mdl_2 = fitlm(cool_temp2,cool_C_2);

cool_mdl_2 = cool_mdl_2.Coefficients;
cool_int_2 = table2array(cool_mdl_2(1,1));
cool_slope_2 = table2array(cool_mdl_2(2,1));

lc2_x = 50:1:72;
lc2_y = cool_int_2 + cool_slope_2*lc2_x;
%
heat_temp2 = ness_temp2(pwm2~=0);
heat_C_2 = C_2(pwm2~=0);
heat_mdl_2 = fitlm(heat_temp2,heat_C_2);

heat_mdl_2 = heat_mdl_2.Coefficients;
heat_int_2 = table2array(heat_mdl_2(1,1));
heat_slope_2 = table2array(heat_mdl_2(2,1));

lh2_x = 50:1:72;
lh2_y = heat_int_2 + heat_slope_2*lh2_x;
%
cool_temp3 = ness_temp3(pwm3==0);
cool_C_3 = C_3(pwm3==0);
cool_mdl_3 = fitlm(cool_temp3,cool_C_3);

cool_mdl_3 = cool_mdl_3.Coefficients;
cool_int_3 = table2array(cool_mdl_3(1,1));
cool_slope_3 = table2array(cool_mdl_3(2,1));

lc3_x = 50:1:112;
lc3_y = cool_int_3 + cool_slope_3*lc3_x;
%
heat_temp3 = ness_temp3(pwm3~=0);
heat_C_3 = C_3(pwm3~=0);
heat_mdl_3 = fitlm(heat_temp3,heat_C_3);


heat_mdl_3 = heat_mdl_3.Coefficients;
heat_int_3 = table2array(heat_mdl_3(1,1));
heat_slope_3 = table2array(heat_mdl_3(2,1));

lh3_x = 50:1:112;
lh3_y = heat_int_3 + heat_slope_3*lh3_x;
%
cool_temp4 = ness_temp4(pwm4==0);
cool_C_4 = C_4(pwm4==0);
cool_mdl_4 = fitlm(cool_temp4,cool_C_4);

cool_mdl_4 = cool_mdl_4.Coefficients;
cool_int_4 = table2array(cool_mdl_4(1,1));
cool_slope_4 = table2array(cool_mdl_4(2,1));

lc4_x = 50:1:105;
lc4_y = cool_int_4 + cool_slope_4*lc4_x;
%
heat_temp4 = ness_temp4(pwm4~=0);
heat_C_4 = C_4(pwm4~=0);
heat_mdl_4 = fitlm(heat_temp4,heat_C_4);

heat_mdl_4 = heat_mdl_4.Coefficients;
heat_int_4 = table2array(heat_mdl_4(1,1));
heat_slope_4 = table2array(heat_mdl_4(2,1));

lh4_x = 50:1:105;
lh4_y = heat_int_4 + heat_slope_4*lh4_x;
%
global cool_int_map;
cool_int_map = [cool_int_1,cool_int_2,cool_int_3,cool_int_4];

global cool_slope_map;
cool_slope_map = [cool_slope_1,cool_slope_2,cool_slope_3,cool_slope_3];

global heat_int_map;
heat_int_map = [heat_int_1,heat_int_2,heat_int_3,heat_int_4];

global heat_slope_map;
heat_slope_map = [heat_slope_1,heat_slope_2,heat_slope_3,heat_slope_3];

%% Plotting
picturewidth_singlecolumn = 20; % set this parameter and keep it forever
picturewidth_doublecolumn = 40;
font_size = 20;%17
line_width = 3;
% single_column = true;
% hw_ratio = 0.9; % feel free to play with this ratio 
single_column = false;
hw_ratio = 0.2; % feel free to play with this ratio 
fname = 'C_calibration_k_matrix';
hfig = figure(1);
color_c = [0.3010 0.7450 0.9330]; % sky blue
color_h = [0.9290 0.6940 0.1250]; % yellow/orange
color_purple = [0.4940 0.1840 0.5560];
color_darkblue = [0 0.4470 0.7410];
color_green = [0.4660 0.6740 0.1880];
subplot(1,4,1)
plot(cool_temp1,cool_C_1,'.b',MarkerSize=20)
hold on
plot(heat_temp1,heat_C_1,'.r',MarkerSize=20)
% plot quadratic estimation of c
plot(lc1_x,lc1_y,'-',LineWidth=3,Color=color_c)
plot(lh1_x,lh1_y,'-',LineWidth=3,Color=color_h)
max_temp = max(lc1_x); min_temp = min(lc1_x);
mid_point_temp = min_temp + (max_temp-min_temp)/2;
c_constant_cool = cool_int_1 + cool_slope_1*mid_point_temp;
c_constant_heat = heat_int_1 + heat_slope_1*mid_point_temp;
% plot linear estimation of c
plot(lc1_x,c_constant_cool*ones(length(lc1_x)),'--',LineWidth=3,Color=color_green);
plot(lh1_x,c_constant_heat*ones(length(lh1_x)),'--',LineWidth=3,Color=color_purple);
hold off
% title('Derived Constant for SMA 1')
xlabel('$T\ (^{\circ}C)$')
ylabel('$\alpha_1 \ (N/^{\circ}C)$')

subplot(1,4,2)
plot(cool_temp2,cool_C_2,'.b',MarkerSize=20)
hold on
plot(heat_temp2,heat_C_2,'.r',MarkerSize=20)
% plot quadratic estimation of c
plot(lc2_x,lc2_y,'-',LineWidth=3,Color=color_c)
plot(lh2_x,lh2_y,'-',LineWidth=3,Color=color_h)
max_temp = max(lc2_x); min_temp = min(lc2_x);
mid_point_temp = min_temp + (max_temp-min_temp)/2;
c_constant_cool = cool_int_2 + cool_slope_2*mid_point_temp;
c_constant_heat = heat_int_2 + heat_slope_2*mid_point_temp;
% plot linear estimation of c
plot(lc2_x,c_constant_cool*ones(length(lc2_x)),'--',LineWidth=3,Color=color_green);
plot(lh2_x,c_constant_heat*ones(length(lh2_x)),'--',LineWidth=3,Color=color_purple);
hold off
% title('Derived Constant for SMA 2')
xlabel('$T \ (^{\circ}C)$')
ylabel('$\alpha_2 \ (N/^{\circ}C)$')

subplot(1,4,3)
plot(cool_temp3,cool_C_3,'.b',MarkerSize=20)
hold on
plot(heat_temp3,heat_C_3,'.r',MarkerSize=20)
% plot quadratic estimation of c
plot(lc3_x,lc3_y,'-',LineWidth=3,Color=color_c)
plot(lh3_x,lh3_y,'-',LineWidth=3,Color=color_h)
max_temp = max(lc3_x); min_temp = min(lc3_x);
mid_point_temp = min_temp + (max_temp-min_temp)/2;
c_constant_cool = cool_int_3 + cool_slope_3*mid_point_temp;
c_constant_heat = heat_int_3 + heat_slope_3*mid_point_temp;
% plot linear estimation of c
plot(lc3_x,c_constant_cool*ones(length(lc3_x),1),'--',LineWidth=3,Color=color_green);
plot(lh3_x,c_constant_heat*ones(length(lh3_x),1),'--',LineWidth=3,Color=color_purple);
hold off
% title('Derived Constant for SMA 3')
xlabel('$T \ (^{\circ}C)$')
ylabel('$\alpha_3 \ (N/^{\circ}C)$')

subplot(1,4,4)
plot(cool_temp4,cool_C_4,'.b',MarkerSize=20)
hold on
plot(heat_temp4,heat_C_4,'.r',MarkerSize=20)
% plot quadratic estimation of c
plot(lc4_x,lc4_y,'-',LineWidth=3,Color=color_c)
plot(lh4_x,lh4_y,'-',LineWidth=3,Color=color_h)
max_temp = max(lc4_x); min_temp = min(lc4_x);
mid_point_temp = min_temp + (max_temp-min_temp)/2;
c_constant_cool = cool_int_4 + cool_slope_4*mid_point_temp;
c_constant_heat = heat_int_4 + heat_slope_4*mid_point_temp;
% plot linear estimation of c
plot(lc4_x,c_constant_cool*ones(length(lc4_x),1),'--',LineWidth=3,Color=color_green);
plot(lh4_x,c_constant_heat*ones(length(lh4_x),1),'--',LineWidth=3,Color=color_purple);
hold off
% title('Derived Constant for SMA 4')
xlabel('$T \ (^{\circ}C)$')
ylabel('$\alpha_4 \ (N/^{\circ}C)$')

set(findall(hfig,'-property','FontSize'),'FontSize',font_size) % adjust fontsize to your document
set(findall(hfig,'-property','LineWidth'),'LineWidth',line_width)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

if single_column==true
    set(hfig,'Units','centimeters','Position',[3 3 picturewidth_singlecolumn hw_ratio*picturewidth_singlecolumn]);
else
    set(hfig,'Units','centimeters','Position',[3 3 picturewidth_doublecolumn hw_ratio*picturewidth_doublecolumn]);
end
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig, strcat('./output_figures/',fname,'.pdf'), 'ContentType', 'vector');
exportgraphics(hfig, strcat('./output_figures/',fname,'.png'), 'ContentType', 'vector');

%%
C_constants_heat =  zeros(4,1);
C_constants_cool =  zeros(4,1);

max_temp = max(lc1_x); min_temp = min(lc1_x);
mid_point_temp = min_temp + (max_temp-min_temp)/2;
C_constants_cool(1) = cool_int_1 + cool_slope_1*mid_point_temp;
C_constants_heat(1) = heat_int_1 + heat_slope_1*mid_point_temp;
max_temp = max(lc2_x); min_temp = min(lc2_x);
mid_point_temp = min_temp + (max_temp-min_temp)/2;
C_constants_cool(2) = cool_int_2 + cool_slope_2*mid_point_temp;
C_constants_heat(2) = heat_int_2 + heat_slope_2*mid_point_temp;
max_temp = max(lc3_x); min_temp = min(lc3_x);
mid_point_temp = min_temp + (max_temp-min_temp)/2;
C_constants_cool(3)  = cool_int_3 + cool_slope_3*mid_point_temp;
C_constants_heat(3) = heat_int_3 + heat_slope_3*mid_point_temp;
max_temp = max(lc4_x); min_temp = min(lc4_x);
mid_point_temp = min_temp + (max_temp-min_temp)/2;
C_constants_cool(4)  = cool_int_4 + cool_slope_4*mid_point_temp;
C_constants_heat(4) = heat_int_4 + heat_slope_4*mid_point_temp;

save('./data4plot/Kmatrix/C_constants_heat_k_matrix.mat', 'C_constants_heat');
save('./data4plot/kmatrix/C_constants_cool_k_matrix.mat', 'C_constants_cool');

%% ONLY FOR the LEGEND
fname = "legend_gen";
hfig = figure(2);
color_c = [0.3010 0.7450 0.9330]; % sky blue
color_h = [0.9290 0.6940 0.1250]; % yellow/orange
color_purple = [0.4940 0.1840 0.5560];
color_darkblue = [0 0.4470 0.7410];
color_green = [0.4660 0.6740 0.1880];

plot(cool_temp1,cool_C_1,'.b',MarkerSize=20)
hold on
plot(heat_temp1,heat_C_1,'.r',MarkerSize=20)
% plot quadratic estimation of c
plot(lc1_x,lc1_y,'-',LineWidth=3,Color=color_c)
plot(lh1_x,lh1_y,'-',LineWidth=3,Color=color_h)
max_temp = max(lc1_x); min_temp = min(lc1_x);
mid_point_temp = min_temp + (max_temp-min_temp)/2;
c_constant_cool = cool_int_1 + cool_slope_1*mid_point_temp;
c_constant_heat = heat_int_1 + heat_slope_1*mid_point_temp;
% plot linear estimation of c
plot(lc1_x,c_constant_cool*ones(length(lc1_x),1),'--',LineWidth=3,Color=color_green);
plot(lh1_x,c_constant_heat*ones(length(lh1_x),1),'--',LineWidth=3,Color=color_purple);
hold off
% title('Derived Constant for SMA 1')
xlabel('$T\ (^{\circ}C)$')
ylabel('$C_1 \ (N/^{\circ}C)$')

legends = {'$Cooling$','$Heating$','$C_{BF}(T)$','$H_{BF}(T)$','$H_{MP}$','$C_{MP}$'};
leg = legend(legends);
legend('Location','eastoutside'); 
leg.ItemTokenSize = [40,40];

set(findall(hfig,'-property','FontSize'),'FontSize',font_size) % adjust fontsize to your document
set(findall(hfig,'-property','LineWidth'),'LineWidth',line_width)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(leg, 'Interpreter','latex');
if single_column==true
    set(hfig,'Units','centimeters','Position',[3 3 picturewidth_singlecolumn hw_ratio*picturewidth_singlecolumn]);
else
    set(hfig,'Units','centimeters','Position',[3 3 picturewidth_doublecolumn hw_ratio*picturewidth_doublecolumn]);
end
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig, strcat('./output_figures/',fname,'.pdf'), 'ContentType', 'vector');
exportgraphics(hfig, strcat('./output_figures/',fname,'.png'), 'ContentType', 'vector');