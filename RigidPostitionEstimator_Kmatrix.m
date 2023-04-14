
% AOI is the angle of interest -- to fully for ne wangles, you need the
% starting angles of the arm when no force is applied

num_links = 8;
AOI = 4;

total_length = 0.188; % total length of limb in meters
total_mass = 0.115; % total mass of limb in kilograms
limb_mass = total_mass/num_links;% mass of each rigid link
limb_length = total_length/num_links; % length of each rigid link 
Length = ones(1,num_links)*limb_length; % input length vector for Jacobian --- vector in case different lengths for links is wanted 

filedir = "./Test-data/test_data_04_04"; % file name of data 
opts = detectImportOptions(filedir,'NumHeaderLines',3); % number of header lines which are to be ignored
opts.VariableNamesLine = 3; % row number which has variable names
% opts.ReadVariableNames = true;
% opts.PreserveVariableNames = true;
% opts.DataLine = 7; % row number from which the actual data starts
% opts.VariableUnitsLine = 6; % row number which has units specified
data = readtable(filedir,opts);
% data = data_original(2:651,:); %% only keep the first few rows from second 0 to 130
%%
data = data(8:end,:);
name_list = data.Properties.VariableNames;

time = table2array(data(:,"TestTime"));
% setpoint = data(:,contains(name_list,'setpoint'));
% setpoint_array = table2array(setpoint);
TOI = find(abs(time - AOI) <= 0.1); % finds point of interest based on angle of interest 
TOI = TOI(1);

% pull out the correct columns from the large table
%PWM = data(:,contains(name_list,'pwm'));
%temp_limit = data(:,contains(name_list,'temp_limit'));
temp = data(TOI,contains(name_list,'MCPtemp'));
temp_start = data(1,contains(name_list,'MCPtemp'));
%angle = data(TOI,contains(name_list,'Bendlabs'));
%setpoint = data(:,contains(name_list,'setpoint'));
tag_position = data(TOI,contains(name_list,'tag'));

all_tag_x = data(:,contains(name_list,'tagid_x'));
all_tag_x_tab = table2array(all_tag_x);
all_tag_y = data(:,contains(name_list,'tagid_y'));
all_tag_y_tab = table2array(all_tag_y);
for i = 1 : length(all_tag_x_tab)
    all_tag_x_tab(i,:) = all_tag_x_tab(i,:) - all_tag_x_tab(i,1);
    all_tag_y_tab(i,:) = all_tag_y_tab(i,:) - all_tag_y_tab(i,1);
end 
end_eff = [(all_tag_x_tab(:,end)-all_tag_x_tab(:,1))' ; (all_tag_y_tab(:,end)-all_tag_y_tab(:,1))'];

tag_position_start = data(1,contains(name_list,'tag'));
pwm_data = data(TOI,contains(name_list,'pwm_duty_limb'));
pwm_data = table2array(pwm_data);

% convert from table to arrays 
%PWM_array = table2array(PWM);
%temp_limit = table2array(temp_limit);
temp = table2array(temp);
initial_temp = table2array(temp_start);
%angle_array = table2array(angle);
%setpoint_array = table2array(setpoint);
tag_position_array = table2array(tag_position);
tag_array(1,:) = tag_position_array(1:2:end);
tag_array(2,:) = -1*tag_position_array(2:2:end);
tag_position_array_start =table2array(tag_position_start);
tag_array_start(1,:) = tag_position_array_start(1:2:end);
tag_array_start(2,:) = tag_position_array_start(2:2:end);

%% conversion variable for pixel to meters
a = zeros(1,num_links);
b = zeros(1,num_links);
convert = zeros(1,num_links);
convert_array_x = tag_array_start(1,:) - tag_array_start(1,1);
convert_array_x = [convert_array_x(1) convert_array_x((8/num_links)+1:(8/num_links):end)];
convert_array_y = tag_array_start(2,:) - tag_array_start(2,1);
convert_array_y = [convert_array_y(1) convert_array_y((8/num_links)+1:(8/num_links):end)];
convert_mat = [convert_array_x ; convert_array_y];

for ii = 1 : num_links

    a(ii) = convert_array_x(ii+1)-convert_array_x(ii);
    
    b(ii) = convert_array_y(ii+1)-convert_array_y(ii);
    
    convert(ii) = (total_length)/(num_links*sqrt(a(ii)^2+b(ii)^2));

end 

convert = mean(convert);

%% Pixel position in meters for comparison

PixPos = PixelToPosition(convert_mat,convert);
PixPos = [PixPos(1,:) ; -1*PixPos(2,:)];
% PixPos = [PixPos(1,:)-PixPos(1,1) ; PixPos(2,:)-PixPos(2,1)];
x = PixPos(1,:);
y = PixPos(2,:);

% plotting parameter is here
single_column = true;
picturewidth_singlecolumn = 20; % set this parameter and keep it forever
picturewidth_doublecolumn = 40;
hw_ratio = 1; % feel free to play with this ratio 
font_size = 20;%17
line_width = 3;
fname = 'Pixel_to_position_8';
% legends = {'Leg1','Leg2','Leg3','Spine1','Spine2','Setpoints'};
hfig= figure(1);  % save the figure handle in a variable
hold on;

% plotting ... customized here
plot(x, y,'.k-',LineWidth=3,MarkerSize=30);
axis equal
axis([-0.01 0.2 -0.11 0.01])
hold on

%
xlabel('$x (m)$')
ylabel('$y (m)$')

% set(gca,'XTick',[0:25:125]);
% set(gca,'xticklabel',[])
% set(gca,'yTick',[-50:25:50]);

% leg = legend(legends);
% legend('Location','eastoutside');     
% leg.ItemTokenSize = [40,40];
% leg.FontSize = font_size;

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
exportgraphics(hfig, strcat(fname,'.pdf'), 'ContentType', 'vector');
exportgraphics(hfig, strcat(fname,'.png'), 'ContentType', 'vector');

% saveas(figure(1),'PixelToPositionPlot_8.png')

% calculates angles between april tags
calc_angles = angle_solver(tag_array);
calc_angles_start = angle_solver(tag_array_start);

% calculating external forces --- only gravity in this case 
g = 9.8; % gravity
alpha_e = [0;limb_mass*g]; 

%% solving for spring constant
% initializing the analytical jacobian 
jacobian_a = zeros(num_links,1);
num = length(calc_angles)/num_links;
angles_of_interest = calc_angles(1:num:end);
angles_of_interest_start = calc_angles_start(1:num:end);

% [x,y] = pose_plotter(-angles_of_interest_start,Length);
% figure()
% plot(x,y,'-sc',LineWidth=2)
    
    for iii = 1 : num_links
        
        link_jacobian_a = transpose(LinkIndependentFormulation(iii,num_links,Length,angles_of_interest_start))*alpha_e;
        jacobian_a = jacobian_a + link_jacobian_a;

    end 

jacobian_a = double(jacobian_a);

% set rest angles 
rest_angles = angles_of_interest_start; 
rest_angles_deg = rad2deg(rest_angles);

% display time step of interest and starting angles
disp('Time in seconds : ')
disp(num2str(time(TOI))); % not necessary for final 
disp('Starting angles in degrees: ')
disp(num2str(rest_angles_deg)); % should be constant

% solve for spring constant
% eqn = jacobian_a == k.*()

for i = 1 : length(jacobian_a)
    k_tor_spring_rest(i) = jacobian_a(i)/rest_angles(i);
end 
% k_tor_spring_rest = linsolve(-jacobian_a,rest_angles');

k_tor_spring_rest = double(k_tor_spring_rest);
% no_grav_angles = zeros(1,length(angles_of_interest_start));
% k_tor_spring_rest = linsolve(-jacobian_a,angles_of_interest_start' - no_grav_angles');
% k_tor_spring_rest = pinv(angles_of_interest' - angles_of_interest_start')*(jacobian_a); % all angles must be same units 
% k_tor_spring_rest = double(k_tor_spring_rest);

disp('Derived Torsional Spring Constant at rest: ')
disp(k_tor_spring_rest)

% solving for pose
pose_angles = pose_solver(k_tor_spring_rest,jacobian_a);

% plotting the derived pose 

[x,y] = pose_plotter(pose_angles,Length);
kmatrix_grav_est = [x;y];
plot(x,y,'-p',LineWidth=2)

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
exportgraphics(hfig, strcat('./figures/',fname,'.pdf'), 'ContentType', 'vector');
exportgraphics(hfig, strcat('./figures/',fname,'.png'), 'ContentType', 'vector');

% Error between theoretical and actual 
percent_error = zeros(1,length(pose_angles));

    for i = 1 : length(pose_angles) 
        percent_error(i) = abs((calc_angles_start(i) - (-1*pose_angles(i)))/(-1*pose_angles(i))) * 100;
    end 

disp('Percent error from derived angles : ')
disp(percent_error)
disp('Initial temperature : ')
disp(initial_temp)
disp('Current temperature : ')
disp(temp)

% Error between to bendlab angles and actual
% percent_error_bendlab = zeros(1,length(pose_angles));
% 
%     for i = 1 : length(pose_angles) 
%         percent_error(i) = abs((calc_angles_start(i) - pose_angles(i))/pose_angles(i)) * 100;
%     end 
% 
% disp('Percent error from bending sensor : ', percent_error)

%% setting up torque matrix for C solving
load('./Previous-results/CaseByCase/Kmatrix/C_1.mat');
load('./Previous-results/CaseByCase/Kmatrix/C_2.mat');
load('./Previous-results/CaseByCase/Kmatrix/C_3.mat');
load('./Previous-results/CaseByCase/Kmatrix/C_4.mat');
load('./Previous-results/ness_temp1.mat');
load('./Previous-results/ness_temp2.mat');
load('./Previous-results/ness_temp3.mat');
load('./Previous-results/ness_temp4.mat');
load('./Previous-results/pwm1.mat');
load('./Previous-results/pwm2.mat');
load('./Previous-results/pwm3.mat');
load('./Previous-results/pwm4.mat');

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


% [C_1,ness_temp1] = independent_C_solver(1,"./SMA-data/test_c1_only.csv",num_links,Length,k_tor_spring_rest,alpha_e);
single_column = false;
hw_ratio = 0.2; % feel free to play with this ratio 
hfig2 = figure(2);
color_c = [0.3010 0.7450 0.9330]; % sky blue
color_h = [0.9290 0.6940 0.1250]; % yellow/orange
subplot(1,4,1)
plot(cool_temp1,cool_C_1,'.b',MarkerSize=20)
hold on
plot(heat_temp1,heat_C_1,'.r',MarkerSize=20)
plot(lc1_x,lc1_y,'-',LineWidth=3,Color=color_c)
plot(lh1_x,lh1_y,'-',LineWidth=3,Color=color_h)
hold off
% title('Derived Constant for SMA 1')
xlabel('$T\ (^{\circ}C)$')
ylabel('$C_1 \ (N/^{\circ}C)$')

% [C_2,ness_temp2] = independent_C_solver(2,"./SMA-data/test_c2_only.csv",num_links,Length,k_tor_spring_rest,alpha_e);
subplot(1,4,2)
plot(cool_temp2,cool_C_2,'.b',MarkerSize=20)
hold on
plot(heat_temp2,heat_C_2,'.r',MarkerSize=20)
plot(lc2_x,lc2_y,'-',LineWidth=3,Color=color_c)
plot(lh2_x,lh2_y,'-',LineWidth=3,Color=color_h)
hold off
% title('Derived Constant for SMA 2')
xlabel('$T \ (^{\circ}C)$')
ylabel('$C_2 \ (N/^{\circ}C)$')

% [C_3,ness_temp3] = independent_C_solver(3,"./SMA-data/test_c3_only.csv",num_links,Length,k_tor_spring_rest,alpha_e);
subplot(1,4,3)
plot(cool_temp3,cool_C_3,'.b',MarkerSize=20)
hold on
plot(heat_temp3,heat_C_3,'.r',MarkerSize=20)
plot(lc3_x,lc3_y,'-',LineWidth=3,Color=color_c)
plot(lh3_x,lh3_y,'-',LineWidth=3,Color=color_h)
hold off
% title('Derived Constant for SMA 3')
xlabel('$T \ (^{\circ}C)$')
ylabel('$C_3 \ (N/^{\circ}C)$')

% [C_4,ness_temp4] = independent_C_solver(4,"./SMA-data/test_c4_only.csv",num_links,Length,k_tor_spring_rest,alpha_e);
subplot(1,4,4)
plot(cool_temp4,cool_C_4,'.b',MarkerSize=20)
hold on
plot(heat_temp4,heat_C_4,'.r',MarkerSize=20)
plot(lc4_x,lc4_y,'-',LineWidth=3,Color=color_c)
plot(lh4_x,lh4_y,'-',LineWidth=3,Color=color_h)
hold off
% title('Derived Constant for SMA 4')
xlabel('$T \ (^{\circ}C)$')
ylabel('$C_4 \ (N/^{\circ}C)$')

fname2 = 'SMA_Linear_Regression';
set(findall(hfig2,'-property','FontSize'),'FontSize',font_size) % adjust fontsize to your document
set(findall(hfig2,'-property','LineWidth'),'LineWidth',line_width)
set(findall(hfig2,'-property','Box'),'Box','off') % optional
set(findall(hfig2,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
if single_column==true
    set(hfig2,'Units','centimeters','Position',[3 3 picturewidth_singlecolumn hw_ratio*picturewidth_singlecolumn]);
else
    set(hfig2,'Units','centimeters','Position',[3 3 picturewidth_doublecolumn hw_ratio*picturewidth_doublecolumn]);
end
pos = get(hfig2,'Position');
set(hfig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig2, strcat('./figures/',fname2,'.pdf'), 'ContentType', 'vector');
exportgraphics(hfig2, strcat('./figures/',fname2,'.png'), 'ContentType', 'vector');

% indx1 = find(C_1>0.28);
% C_1 = C_1(indx1);
% ness_temp1 = ness_temp1(indx1);
% indx2 = find(C_2>0.43);
% C_2 = C_2(indx2);
% ness_temp2 = ness_temp2(indx2);
% indx3 = find(C_3>0.49);
% C_3 = C_3(indx3);
% ness_temp3 = ness_temp3(indx3);
% indx4 = find(C_4>0.5);
% C_4 = C_4(indx4);
% ness_temp4 = ness_temp4(indx4);
% 
% % [C_1,ness_temp1] = independent_C_solver(1,"test_c1_only.csv",num_links,Length,k_tor_spring_rest,alpha_e);
% figure()
% subplot(2,2,1)
% plot(ness_temp1,C_1,'*b')
% title('Derived Constant for SMA 1')
% xlabel('Temperature (Celsius)')
% ylabel('SMA constant')
% % 
% % [C_2,ness_temp2] = independent_C_solver(2,"test_c2_only.csv",num_links,Length,k_tor_spring_rest,alpha_e);
% subplot(2,2,2)
% plot(ness_temp2,C_2,'*b')
% title('Derived Constant for SMA 2')
% xlabel('Temperature (Celsius)')
% ylabel('SMA constant')
% % 
% % [C_3,ness_temp3] = independent_C_solver(3,"test_c3_only.csv",num_links,Length,k_tor_spring_rest,alpha_e);
% subplot(2,2,3)
% plot(ness_temp3,C_3,'*b')
% title('Derived Constant for SMA 3')
% xlabel('Temperature (Celsius)')
% ylabel('SMA constant')
% % 
% % [C_4,ness_temp4] = independent_C_solver(4,"test_c4_only.csv",num_links,Length,k_tor_spring_rest,alpha_e);
% subplot(2,2,4)
% plot(ness_temp4,C_4,'*b')
% title('Derived Constant for SMA 4')
% xlabel('Temperature (Celsius)')
% ylabel('SMA constant')
%%

% all_temp = data(:,contains(name_list,'MCPtemp'));
% all_temp = table2array(all_temp);
% all_angles = zeros(length(all_temp),8);
% 
%     for bruh = 1 : length(all_temp)
%     
%         tag_position = data(bruh,contains(name_list,'tag'));
%         tag_position_array = table2array(tag_position);
%         tag_array(1,:) = tag_position_array(1:2:end);
%         tag_array(2,:) = -1*tag_position_array(2:2:end);
%         all_angles(bruh,:) = angle_solver(tag_array);
%     
%     end
% 
% %% setting up torque caused by SMA
% 
% [C_matrix] = c_solver(all_temp,all_angles,k_tor_spring_rest,num_links,Length,alpha_e);
% 
% valid_C_tuple = cell(1,4);
% valid_temp_tuple = cell(1,4);
% 
%     for i = 1 : 4
% 
%         idxes = find(all_temp(:,i)>40);
%         valid_C = C_matrix(i,idxes);
%         valid_temp = all_temp(idxes,i);
%         valid_C_tuple{i} = valid_C;
%         valid_temp_tuple{i} = valid_temp;
% 
%     end
% 
% [C_matrix] = [mean(valid_C_tuple{1});mean(valid_C_tuple{2});mean(valid_C_tuple{3});mean(valid_C_tuple{4})];
    
%%

jacobian_a2 = zeros(num_links,1);

    for iii = 1 : num_links
        
        link_jacobian_a = transpose(LinkIndependentFormulation(iii,num_links,Length,angles_of_interest))*alpha_e;
        jacobian_a2 = jacobian_a2 + link_jacobian_a;

    end 

jacobian_a2 = double(jacobian_a2);

final_angles = pose_solver_C(k_tor_spring_rest,jacobian_a2,temp,initial_temp,num_links,limb_length,pwm_data);

disp(rad2deg(final_angles))
%%
hfig3 = figure(3);
[x1,y1] = pose_plotter(calc_angles,Length);
plot(x1,y1,'-*k',LineWidth=2)
axis equal
axis([0 0.2 -0.15 0.01])
hold on
[x2,y2] = pose_plotter(final_angles,Length);
plot(x2,y2,'-*g',LineWidth=2)
xlabel('$x (m)$')
ylabel('$y (m)$')


single_column = true;
hw_ratio = 1;
fname3 = 'Example_Frame_Estimation';
set(findall(hfig3,'-property','FontSize'),'FontSize',font_size) % adjust fontsize to your document
set(findall(hfig3,'-property','LineWidth'),'LineWidth',line_width)
set(findall(hfig3,'-property','Box'),'Box','off') % optional
set(findall(hfig3,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig3,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
if single_column==true
    set(hfig3,'Units','centimeters','Position',[3 3 picturewidth_singlecolumn hw_ratio*picturewidth_singlecolumn]);
else
    set(hfig3,'Units','centimeters','Position',[3 3 picturewidth_doublecolumn hw_ratio*picturewidth_doublecolumn]);
end
pos = get(hfig3,'Position');
set(hfig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
exportgraphics(hfig3, strcat('./figures/',fname3,'.pdf'), 'ContentType', 'vector');
exportgraphics(hfig3, strcat('./figures/',fname3,'.png'), 'ContentType', 'vector');
%%
% pwm_data_all = data(:,contains(name_list,'pwm_duty_limb'));
% pwm_data_all = table2array(pwm_data_all);
% temp_all = data(:,contains(name_list,'MCPtemp'));
% temp_all = table2array(temp_all);
% initial_temp = temp_all(1,:);
% all_final_angles = zeros(length(all_tag_x_tab),num_links);
% all_angles_pixel = zeros(length(all_tag_x_tab),num_links);
% FF = waitbar(0,'Sample indx');
% 
% for I = 1 : length(all_tag_x_tab)
%             waitbar(I/length(all_tag_x_tab),FF,strcat('i = ',sprintf('%d',I)))
%     location_list = [all_tag_x_tab(I,:) ; all_tag_y_tab(I,:)];
%     angles = angle_solver(location_list);
%     all_angles_pixel(I,:) = angles;
%     num = length(calc_angles)/num_links;
%     angles_of_interest = angles(1:num:end);
%     jacobian_a2 = zeros(num_links,1);
% 
%     for iii = 1 : num_links
%         
%         link_jacobian_a = transpose(LinkIndependentFormulation(iii,num_links,Length,angles_of_interest))*alpha_e;
%         jacobian_a2 = jacobian_a2 + link_jacobian_a;
% 
%     end 
% 
%     jacobian_a2 = double(jacobian_a2);
%     
%     temp_AP = temp_all(I,:);
%     pwm_data = pwm_data_all(I,:);
% 
%     all_final_angles(I,:) = pose_solver_C(k_tor_spring_rest,jacobian_a2,temp_AP,initial_temp,num_links,limb_length,pwm_data);
% 
% end
% 
% close(FF)

%%
%%
load('./CaseByCase/Kmatrix/pixel_angles_Cvariable.mat');
load('./CaseByCase/Kmatrix/final_angles_Cvariable.mat');
% load('./CaseByCase/Kmatrix/pixel_angles_Cconstant.mat');
% load('./CaseByCase/Kmatrix/final_angles_Cconstant.mat');
error_x = zeros(1,length(all_angles_pixel));
error_y = error_x;
X1 = zeros(length(all_angles_pixel),num_links+1);
Y1 = X1;
X2 = X1;
Y2 = X1;

for i = 1 : length(all_angles_pixel)

    pix_theta = all_angles_pixel(i,:);
    angle_theta = all_final_angles(i,:);

    [x1,y1] = pose_plotter(-pix_theta,Length);
    X1(i,:) = x1;
    Y1(i,:) = y1;
    figure(4)
    pause(0.1)
    plot(x1,y1,'k')
    hold on
    [x2,y2] = pose_plotter(angle_theta,Length);
    X2(i,:) = x2;
    Y2(i,:) = y2;
    plot(x2,y2,'r')
    hold off

    error_x(i) = abs(x1(end)-x2(end));
    error_y(i) = abs(y1(end)-y2(end));

end
error_mat = [error_x ; error_y];
diff_angle = all_final_angles-all_angles_pixel;
mean_angle4step = mean(diff_angle,2);
std_angle4step = std(diff_angle,0,2);


x = 1 : length(all_angles_pixel);
curve1 = (mean_angle4step+std_angle4step)';
curve2 = (mean_angle4step-std_angle4step)';
% plot(x, curve1, 'k', 'LineWidth', 2);
hold on;
% plot(x, curve2, 'k', 'LineWidth', 2);
plot(x,mean_angle4step,'r', 'LineWidth', 2)
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'black','FaceAlpha',0.2);
hold off

% for i = 1 : length(X1)
%     
%     for ii = 1 : num_links+1
% 
%     a = [X1(i,ii) Y1(i,ii)];
%     b = [X2(i,ii) Y2(i,ii)];
%     d(i,ii) = norm(b-a);
% 
%     end 
% end 
% 
% TTime = 1 : length(X1);
% for i = 1 : num_links+1
% plot(TTime,d(:,i))
% hold on 
% end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function angles = angle_solver(location_list)
% calculated angles are in radians

% calculates number of necessary iterations given by location_list 
iter = length(location_list);

% initilizes vector to old derived angles 
angles = zeros(1,iter-1);
    
    for itr = 1:iter-1

%         x1 = location_list(1,itr);
%         x2 = location_list(1,itr+1);
%         y1 = location_list(2,itr);
%         y2 = location_list(2,itr+1);
% 
%         angles(itr) = atan2(x1*y2-y1*x2,x1*x2+y1*y2);
        
        angles(itr) = atan2((location_list(2,itr+1)-location_list(2,itr)),(location_list(1,itr+1)-location_list(1,itr)));

        if itr == 1 

        else

        angles(itr) = angles(itr)-sum(angles(1:itr-1));

        end 

    end 

end 


function final_angles = pose_solver(spring_constants,jacobian_matrix)
% solves for pose angles based on derived spring constant

% syms final_angles [num_links 1] % could probably do it with '\' instead but this way seems to work without much issue 
%final_angles = zeros(num_links,1);

% eqn = jacobian_matrix - spring_constants*(final_angles - transpose(starting_angles)) == 0;
% final_angles = (jacobian_matrix - spring_constants*transpose(starting_angles))/spring_constants;
final_angles = zeros(1,length(jacobian_matrix));
for i = 1 : length(jacobian_matrix)

final_angles(i) = jacobian_matrix(i)/-spring_constants(i);

end 

% [final_angles] = vpasolve(eqn);

end 

function final_angles = pose_solver_C(spring_constants,jacobian_matrix,temp,initial_temp,num_links,limb_length,pwm_data)
% solves for pose angles based on derived spring constant

% syms final_angles [num_links 1] % could probably do it with '\' instead but this way seems to work without much issue 
%final_angles = zeros(num_links,1);

% eqn = jacobian_matrix - spring_constants*(final_angles - transpose(starting_angles)) == 0;
% final_angles = (jacobian_matrix - spring_constants*transpose(starting_angles))/spring_constants;
C = zeros(num_links,1);
C_matrix = zeros(4,1);

load('./CaseByCase/Kmatrix/C_constants_cool_k_matrix.mat');
C_val_cool = C_constants_cool';
load('./CaseByCase/Kmatrix/C_constants_heat_k_matrix.mat');
C_val_heat = C_constants_heat';

mat = [C_val_heat ; C_val_cool];

% C_matrix(1) = get_c_with_temp(pwm_data(2),temp(1),1);
C_matrix(1) = get_c_from_const(pwm_data(2),1,mat);
% C_matrix(2) = get_c_with_temp(pwm_data(1),temp(2),2);
C_matrix(2) = get_c_from_const(pwm_data(1),2,mat);
% C_matrix(3) = get_c_with_temp(pwm_data(4),temp(4),3);
C_matrix(3) = get_c_from_const(pwm_data(4),3,mat);
% C_matrix(4) = get_c_with_temp(pwm_data(3),temp(3),4);
C_matrix(4) = get_c_from_const(pwm_data(3),4,mat);


% be aware that SMA constant C4 goes with temperatures 3 and C3 goes with
% temperatures 4, this is because we wanted both positive directions of
% movement to be in the same direction 

for ii = 1 : num_links

    if ii <= num_links/2

        C(ii) = (2/num_links*ii*limb_length)*(C_matrix(1)*(temp(1)-initial_temp(1)) - C_matrix(2)*(temp(2)-initial_temp(2)));

    else 

        C(ii) = (2/num_links*ii*limb_length)*(C_matrix(3)*(temp(4)-initial_temp(4)) - C_matrix(4)*(temp(3)-initial_temp(3)));

    end 

end 

mat = jacobian_matrix + C;
final_angles = zeros(1,length(mat));

for i = 1 : length(jacobian_matrix)
final_angles(i) = (-1/spring_constants(i))*mat(i);
end

end 


function [x,y] = pose_plotter(Theta,Length)

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

function PixPos = PixelToPosition(tag_array,convert)

PixPos = tag_array*convert;

end 

function [C_matrix] = c_solver(all_temp,all_angles,k_tor_spring_rest,num_links,Length,alpha_e)

C_matrix = zeros(4,length(all_temp));

limb_length = Length(1);

    for bruh2 = 1 : length(all_temp)
    
        temp = all_temp(bruh2,:);
        initial_temp = all_temp(1,:);

        angles_of_interest_start = all_angles(1,:);
        angles_of_interest = all_angles(bruh2,:);
    
        num = length(angles_of_interest_start)/num_links;
        angles_of_interest = angles_of_interest(1:num:end);
        angles_of_interest_start = angles_of_interest_start(1:num:end);
    
        TWS = k_tor_spring_rest.*(angles_of_interest' - angles_of_interest_start'); % All torque on system excluding SMAs 
    
        
            for iii = 1 : num_links
                
                link_jacobian_a = transpose(LinkIndependentFormulation(iii,num_links,Length,angles_of_interest))*alpha_e;
                TWS = TWS - link_jacobian_a;
        
            end 
    
        TWS = double(TWS);
    
        % initializing the analytical jacobian 
    
        T1 = zeros(num_links/2,2);
        
        T2 = zeros(num_links/2,2);
        
            for itera = 1 : num_links/2
            
                T1(itera,1) = (temp(1)-initial_temp(1))*(2/num_links)*limb_length*itera;
                T1(itera,2) = -1*(temp(2)-initial_temp(2))*(2/num_links)*limb_length*itera;
            
                T2(itera,1) = (temp(4)-initial_temp(4))*(2/num_links)*limb_length*itera;
                T2(itera,2) = -1*(temp(3)-initial_temp(3))*(2/num_links)*limb_length*itera;
            
            end 
        
        %% solving for SMA constants
        
        c1 = pinv(T1)*TWS(1:num_links/2);
        c2 = pinv(T2)*TWS(num_links/2+1:end);
        
        C_matrix(:,bruh2) = [c1;c2];
    
    end 



end 

function [C,ness_temp] = independent_C_solver(C_num,fileid,num_links,Length,k_tor_spring_rest,alpha_e)

opts = detectImportOptions(fileid,'NumHeaderLines',3); % number of header lines which are to be ignored
opts.VariableNamesLine = 3; % row number which has variable names

data_ind = readtable(fileid,opts);

data_ind = data_ind(10:end,:);
name_list_ind = data_ind.Properties.VariableNames;

all_temp_ind = data_ind(:,contains(name_list_ind,'MCPtemp'));
all_temp_ind = table2array(all_temp_ind);
all_angles_ind = zeros(length(all_temp_ind),8);

C = zeros(1,length(all_temp_ind));

limb_length = Length(1);

    for bruh = 1 : length(all_temp_ind)
    
        tag_position_ind = data_ind(bruh,contains(name_list_ind,'tag'));
        tag_position_array_ind = table2array(tag_position_ind);
        tag_array_ind(1,:) = tag_position_array_ind(1:2:end);
        tag_array_ind(2,:) = -1*tag_position_array_ind(2:2:end);
        all_angles_ind(bruh,:) = angle_solver(tag_array_ind);
    
    end
FF = waitbar(0,'Sample indx');

    for bruh2 = 1 : length(all_temp_ind)
        waitbar(bruh2/length(all_temp_ind),FF,strcat('i = ',sprintf('%d',bruh2)))
        temp_ind = all_temp_ind(bruh2,:);
        initial_temp_ind = all_temp_ind(1,:);

        angles_of_interest_start_ind = all_angles_ind(1,:);
        angles_of_interest_ind = all_angles_ind(bruh2,:);
    
        num = length(angles_of_interest_start_ind)/num_links;
        angles_of_interest_ind = angles_of_interest_ind(1:num:end);
        angles_of_interest_start_ind = angles_of_interest_start_ind(1:num:end);
        TWS = zeros(length(angles_of_interest_start_ind),1);

        for i = 1 : length(angles_of_interest_start_ind)

            TWS(i) = k_tor_spring_rest(i)*(angles_of_interest_ind(i)' - angles_of_interest_start_ind(i)'); % All torque on system excluding SMAs 

        end 
    
%         TWS = k_tor_spring_rest.*(angles_of_interest_ind' - angles_of_interest_start_ind'); % All torque on system excluding SMAs 
    
        
            for iii = 1 : num_links
                
                link_jacobian_a = transpose(LinkIndependentFormulation(iii,num_links,Length,angles_of_interest_ind))*alpha_e;
                TWS = TWS + link_jacobian_a;
        
            end 
    
        TWS = double(TWS);
    
        % initializing the analytical jacobian 
    
        T1 = zeros(num_links/2,1);
        
        T2 = zeros(num_links/2,1);

            if C_num == 1
                    
                for itera = 1 : num_links/2
                    T1(itera,1) = (temp_ind(1)-initial_temp_ind(1))*(2/num_links)*limb_length*itera;
                    T1(itera,2) = 0;
                end 

                work_mat = pinv(T1)*-TWS(1:num_links/2);
                c = work_mat(1);
    
            elseif C_num == 2
    
                for itera = 1 : num_links/2

                    T1(itera,1) = 0;
%                     T1(itera,2) = -1*(temp_ind(2)-initial_temp_ind(2))*(2/num_links)*limb_length*itera;
                    T1(itera,2) = (temp_ind(2)-initial_temp_ind(2))*(2/num_links)*limb_length*itera;

                end 

                work_mat = pinv(T1)*-TWS(1:num_links/2);
                c = work_mat(2);
    
            elseif C_num == 3
    
                for itera = 1 : num_links/2

                    T2(itera,1) = (temp_ind(4)-initial_temp_ind(4))*(2/num_links)*limb_length*itera;
                    T2(itera,2) = 0;

                end

                work_mat = pinv(T2)*-TWS(num_links/2+1:end);
                c = work_mat(1);
    
            else
    
                for itera = 1 : num_links/2

                    T2(itera,1) = 0;   
%                     T2(itera,2) = -1*(temp_ind(3)-initial_temp_ind(3))*(2/num_links)*limb_length*itera;
                    T2(itera,2) = (temp_ind(3)-initial_temp_ind(3))*(2/num_links)*limb_length*itera;

                end

                work_mat = pinv(T2)*-TWS(num_links/2+1:end);
                c = work_mat(2);
    
            end 
        
        C(:,bruh2) = c;

    end
    close(FF)

    if C_num == 1
    ness_temp = all_temp_ind(:,C_num);
    elseif C_num == 2
    ness_temp = all_temp_ind(:,C_num);
    elseif C_num == 3 
    ness_temp = all_temp_ind(:,C_num+1);
    elseif C_num == 4
    ness_temp = all_temp_ind(:,C_num-1);
    else
    ness_temp = all_temp_ind(:,C_num);
    end

end 



function  [C] = get_c_with_temp(pwm,temp,sma_num)

    global cool_int_map ;
    global cool_slope_map ;
    global heat_int_map ;
    global heat_slope_map ;

    if (pwm>0)

        C = heat_int_map(sma_num) + temp*heat_slope_map(sma_num);

    else

        C = cool_int_map(sma_num) + temp*cool_slope_map(sma_num);
    
    end

end

function [C] = get_c_from_const(pwm,sma_num,mat)

    if (pwm>0)

        C = mat(1,sma_num);

    else

        C = mat(2,sma_num);
    
    end

end 
