clc;close all;clear;
load('./data4plot/Kconstant/pixel_angles_Cconstant.mat');
all_ground_truth_angles = all_angles_pixel;

load('./data4plot/Kconstant/final_angles_Cconstant.mat');
constant_k_constant_c_all_final_angles = all_final_angles;
load('./data4plot/Kconstant/final_angles_Cvariable.mat');
constant_k_variable_c_all_final_angles = all_final_angles;
load('./data4plot/Kmatrix/final_angles_Cconstant.mat');
k_matrix_constant_c_all_final_angles = all_final_angles;
load('./data4plot/Kmatrix/final_angles_Cvariable.mat');
k_matrix_variable_c_all_final_angles = all_final_angles;

load('./data4plot/Kconstant/time.mat');
time = time';
%%
critical_time = 19.6;
critical_time_idx = find(abs(time - critical_time)<0.1);

num_links = 8;
total_length = 0.188; % total length of limb in meters
limb_length = total_length/num_links; % length of each rigid link 
Length = ones(1,num_links)*limb_length; % input length vector for Jacobian --- vector in case different lengths for links is wanted 
[xf,yf] = pose_plotter_pixel_angle(all_ground_truth_angles(critical_time_idx,:),Length);
[x1,y1] = pose_plotter(constant_k_constant_c_all_final_angles(critical_time_idx,:),Length);
[x2,y2] = pose_plotter(constant_k_variable_c_all_final_angles(critical_time_idx,:),Length);
[x3,y3] = pose_plotter(k_matrix_constant_c_all_final_angles(critical_time_idx,:),Length);
[x4,y4] = pose_plotter(k_matrix_variable_c_all_final_angles(critical_time_idx,:),Length);
%% Plotting
picturewidth_singlecolumn = 20; % set this parameter and keep it forever
picturewidth_doublecolumn = 40;
font_size = 20;%17
line_width = 3;
% single_column = true;
% hw_ratio = 0.9; % feel free to play with this ratio 
single_column = true;
hw_ratio = 1; % feel free to play with this ratio 
fname = 'pose_camparison';
hfig = figure(1);

color_map = parula(8);
color_map = [color_map(1,:);color_map(3,:);color_map(5,:);color_map(7,:)];
hold on
plot(xf, yf,'d-k',MarkerSize=10);
plot(x1 ,y1,'.-','color',color_map(1,:),MarkerSize=30);
plot(x2 ,y2,'.-','color',color_map(2,:),MarkerSize=30);
plot(x3 ,y3,'.-','color',color_map(3,:),MarkerSize=30);
plot(x4 ,y4,'.-','color',color_map(4,:),MarkerSize=30);
hold off
% title('Derived Constant for SMA 1')
xlabel('$x\ (m)$');
ylabel('$y\ (m)$');
axis equal
axis([0 0.2 -0.12 0.06])

legends = {'$CV_{ref}$', '$Case\ 1$','$Case\ 2$','$Case\ 3$','$Case\ 4$'};
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
%% 
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig, strcat('./output_figures/',fname,'.pdf'), 'ContentType', 'vector');
exportgraphics(hfig, strcat('./output_figures/',fname,'.png'), 'ContentType', 'vector');
%%
function [x,y] = pose_plotter(Theta,Length)
Theta = Theta;
xxx = zeros(1,length(Theta));
yyy = zeros(1,length(Theta));
    
    for ii = 1 : length(Theta)

            if ii == 1 
                xxx(ii) = Length(ii)*cos(Theta(ii));
                yyy(ii) = Length(ii)*sin(Theta(ii));
            else 
                xxx(ii) = xxx(ii-1) + Length(ii)*cos(sum(Theta(1:ii)));
                yyy(ii) = yyy(ii-1) + Length(ii)*sin(sum(Theta(1:ii)));
            end 
    
    end 

x = [0 xxx];
y = [0 yyy];

end 
%%
function [x,y] = pose_plotter_pixel_angle(Theta,Length)
Theta = -Theta;
xxx = zeros(1,length(Theta));
yyy = zeros(1,length(Theta));
    
    for ii = 1 : length(Theta)

            if ii == 1 
                xxx(ii) = Length(ii)*cos(Theta(ii));
                yyy(ii) = Length(ii)*sin(Theta(ii));
            else 
                xxx(ii) = xxx(ii-1) + Length(ii)*cos(sum(Theta(1:ii)));
                yyy(ii) = yyy(ii-1) + Length(ii)*sin(sum(Theta(1:ii)));
            end 
    
    end 

x = [0 xxx];
y = [0 yyy];

end 