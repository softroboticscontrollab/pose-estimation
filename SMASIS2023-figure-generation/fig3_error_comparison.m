clc;close all;clear;

load('./data4plot/Kconstant/error_matrix_constant.mat');
constant_k_constant_c_error = error_mat;
load('./data4plot/Kconstant/error_matrix_variable.mat');
constant_k_variable_c_error = error_mat;
load('./data4plot/Kmatrix/error_matrix_constant.mat');
k_matrix_constant_c_error = error_mat;
load('./data4plot/Kmatrix/error_matrix_variable.mat');
k_matrix_variable_c_error = error_mat;
load('./data4plot/Kconstant/time.mat');
time = time';
%%
constant_k_constant_c_error_norm = err2norm(constant_k_constant_c_error);
constant_k_variable_c_error_norm = err2norm(constant_k_variable_c_error);
k_matrix_constant_c_error_norm = err2norm(k_matrix_constant_c_error);
k_matrix_variable_c_error_norm = err2norm(k_matrix_variable_c_error);
%%
%% Plotting
picturewidth_singlecolumn = 20; % set this parameter and keep it forever
picturewidth_doublecolumn = 40;
font_size = 20;%17
line_width = 3;
% single_column = true;
% hw_ratio = 0.9; % feel free to play with this ratio 
single_column = false;
hw_ratio = 0.2; % feel free to play with this ratio 
fname = 'error_camparison';
hfig = figure(1);
color_c = [0.3010 0.7450 0.9330]; % sky blue
color_h = [0.9290 0.6940 0.1250]; % yellow/orange
color_purple = [0.4940 0.1840 0.5560];
color_darkblue = [0 0.4470 0.7410];
color_green = [0.4660 0.6740 0.1880];

color_map = parula(8);
color_map = [color_map(1,:);color_map(3,:);color_map(5,:);color_map(7,:)];
hold on
plot(time ,constant_k_constant_c_error_norm,'-','color',color_map(1,:));
plot(time ,constant_k_variable_c_error_norm,'-','color',color_map(2,:));
plot(time ,k_matrix_constant_c_error_norm,'-','color',color_map(3,:));
plot(time ,k_matrix_variable_c_error_norm,'-','color',color_map(4,:));

critical_time = 19.6;
plot([critical_time,critical_time],[0 0.3],'-r')

hold off
% title('Derived Constant for SMA 1')
xlabel('$t\ (s)$');
ylabel('$||e_{i}||\ (m)$');

legends = {'$Case\ 1$','$Case\ 2$','$Case\ 3$','$Case\ 4$'};
leg = legend(legends);
legend('Location','eastoutside'); 
leg.ItemTokenSize = [40,40];

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
function norm = err2norm(err_matrix)
    [~,sample_size] = size(err_matrix);
    norm = zeros(1,sample_size);
    for idx = 1:sample_size
        norm(idx) = sqrt(err_matrix(1,idx)^2+err_matrix(2,idx)^2);
    end
end
