clc;close all;clear;
num_links = 8;
AOI = 4;

total_length = 0.188; % total length of limb in meters
total_mass = 0.115; % total mass of limb in kilograms
limb_mass = total_mass/num_links;% mass of each rigid link
limb_length = total_length/num_links; % length of each rigid link 
Length = ones(1,num_links)*limb_length; % input length vector for Jacobian --- vector in case different lengths for links is wanted 

filedir = "test_data_04_04"; % file name of data 
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
% a = zeros(1,num_links);
% b = zeros(1,num_links);
% convert = zeros(1,num_links);
% convert_array_x = tag_array_start(1,:) - tag_array_start(1,1);
% convert_array_x = [convert_array_x(1) convert_array_x((8/num_links)+1:(8/num_links):end)];
% convert_array_y = tag_array_start(2,:) - tag_array_start(2,1);
% convert_array_y = [convert_array_y(1) convert_array_y((8/num_links)+1:(8/num_links):end)];
% convert_mat = [convert_array_x ; convert_array_y];
% 
% for ii = 1 : num_links
% 
%     a(ii) = convert_array_x(ii+1)-convert_array_x(ii);
%     
%     b(ii) = convert_array_y(ii+1)-convert_array_y(ii);
%     
%     convert(ii) = (total_length)/(num_links*sqrt(a(ii)^2+b(ii)^2));
% 
% end 
% 
% convert = mean(convert);

%% Pixel position in meters for comparison
% 
% PixPos = PixelToPosition(convert_mat,convert);
% PixPos = [PixPos(1,:) ; -1*PixPos(2,:)];
% % PixPos = [PixPos(1,:)-PixPos(1,1) ; PixPos(2,:)-PixPos(2,1)];
% x = PixPos(1,:);
% y = PixPos(2,:);

% % plotting parameter is here
% single_column = true;
% picturewidth_singlecolumn = 20; % set this parameter and keep it forever
% picturewidth_doublecolumn = 40;
% hw_ratio = 1; % feel free to play with this ratio 
% font_size = 20;%17
% line_width = 3;
% fname = 'Pixel_to_position_8';
% % legends = {'Leg1','Leg2','Leg3','Spine1','Spine2','Setpoints'};
% hfig= figure(1);  % save the figure handle in a variable
% hold on;
% 
% % plotting ... customized here
% plot(x, y,'.k-',LineWidth=3,MarkerSize=30);
% axis equal
% axis([-0.01 0.2 -0.11 0.01])
% hold on
% 
% %
% xlabel('$x (m)$')
% ylabel('$y (m)$')
% 
% % set(gca,'XTick',[0:25:125]);
% % set(gca,'xticklabel',[])
% % set(gca,'yTick',[-50:25:50]);
% 
% % leg = legend(legends);
% % legend('Location','eastoutside');     
% % leg.ItemTokenSize = [40,40];
% % leg.FontSize = font_size;
% 
% set(findall(hfig,'-property','FontSize'),'FontSize',font_size) % adjust fontsize to your document
% set(findall(hfig,'-property','LineWidth'),'LineWidth',line_width)
% set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% if single_column==true
%     set(hfig,'Units','centimeters','Position',[3 3 picturewidth_singlecolumn hw_ratio*picturewidth_singlecolumn]);
% else
%     set(hfig,'Units','centimeters','Position',[3 3 picturewidth_doublecolumn hw_ratio*picturewidth_doublecolumn]);
% end
% pos = get(hfig,'Position');
% set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% exportgraphics(hfig, strcat(fname,'.pdf'), 'ContentType', 'vector');
% exportgraphics(hfig, strcat(fname,'.png'), 'ContentType', 'vector');

%%
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
    k_matrix_tor_spring_rest(i) = jacobian_a(i)/rest_angles(i);
end 
k_matrix_tor_spring_rest = double(k_matrix_tor_spring_rest);
disp('Derived Torsional Spring Constants at rest: ')
disp(k_matrix_tor_spring_rest)

k_constant_tor_spring_rest = linsolve(-jacobian_a,angles_of_interest_start');
k_constant_tor_spring_rest = double(k_constant_tor_spring_rest);

% solving for pose
pose_angles_k_matrix = pose_solver_k_matrix(k_matrix_tor_spring_rest,jacobian_a);
pose_angles_k_constant = pose_solver_k_constant(k_constant_tor_spring_rest,jacobian_a);
%% Convert pixel to meter
%conversion variable for pixel to meters
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

%Pixel position in meters for comparison
PixPos = PixelToPosition(convert_mat,convert);
PixPos = [PixPos(1,:) ; -1*PixPos(2,:)];
% PixPos = [PixPos(1,:)-PixPos(1,1) ; PixPos(2,:)-PixPos(2,1)];
xf = PixPos(1,:);
yf = PixPos(2,:);
%% PLOTTING

single_column = true;
picturewidth_singlecolumn = 20; % set this parameter and keep it forever
picturewidth_doublecolumn = 40;
hw_ratio = 1; % feel free to play with this ratio 
font_size = 20;%17
line_width = 3;
fname = 'rest_pose_comparison';
% legends = {'Leg1','Leg2','Leg3','Spine1','Spine2','Setpoints'};
color_map = parula(8);
color_map = [color_map(1,:);color_map(3,:);color_map(5,:);color_map(7,:)];
hold on
hfig= figure(1);  % save the figure handle in a variable
%[xf,yf] = pose_plotter_pixel_angle(angles_of_interest_start,Length);
[x1,y1] = pose_plotter(pose_angles_k_constant,Length);
[x2,y2] = pose_plotter(pose_angles_k_matrix,Length);
plot(xf, yf,'d-k',MarkerSize=10);
plot(x1 ,y1,'.-','color',color_map(1,:),MarkerSize=30);
plot(x2 ,y2,'.-','color',color_map(2,:),MarkerSize=30);

xlabel('$x\ (m)$');
ylabel('$y\ (m)$');
axis equal
axis([0 0.2 -0.12 0.06])

legends = {'$CV_{ref}$', '$Constant \ K$','$Matrix \ K$'};
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
function angles = angle_solver(location_list)
% calculated angles are in radians

% calculates number of necessary iterations given by location_list 
iter = length(location_list);

% initilizes vector to old derived angles 
angles = zeros(1,iter-1);
    
    for itr = 1:iter-1
        
        angles(itr) = atan2((location_list(2,itr+1)-location_list(2,itr)),(location_list(1,itr+1)-location_list(1,itr)));

        if itr == 1 

        else

        angles(itr) = angles(itr)-sum(angles(1:itr-1));

        end 

    end 

end 


function final_angles = pose_solver_k_matrix(spring_constants,jacobian_matrix)
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

function final_angles = pose_solver_k_constant(spring_constant,jacobian_matrix)
% solves for pose angles based on derived spring constant

% syms final_angles [num_links 1] % could probably do it with '\' instead but this way seems to work without much issue 
% final_angles = zeros(num_links,1);

% eqn = jacobian_matrix - spring_constant*(final_angles - transpose(starting_angles)) == 0;
% final_angles = (jacobian_matrix - spring_constant*transpose(starting_angles))/spring_constant;
% 
% [final_angles] = vpasolve(eqn);
final_angles = jacobian_matrix/-spring_constant;

end 
%%
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