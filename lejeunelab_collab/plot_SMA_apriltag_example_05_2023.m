% Example script that plots the SMA-powered robot arm position, based on
% the AprilTag fiducials, at an example time of interest.
% (C) Soft Robotics Control Lab 2023

clc;close all;clear;

%% Setup and read in data

num_links = 8; % there are eight segments between the fiducial tags. Needed to parse the columns of the .csv file.
total_length = 0.188; % total length of limb in meters. Used to approximately plot the x,y fiducial locations in meters rather than pixel coordinates.

filename = "ezloophw_openloop_datarx_cvtracking_2_limbs_2023-4-5_105056.csv"; % file name of data 
opts = detectImportOptions(filename,'NumHeaderLines',3); % number of header lines which are to be ignored
opts.VariableNamesLine = 3; % row number which has variable names
data = readtable(filename,opts);
data = data(8:end,:); % skip the first handful of rows, junk
name_list = data.Properties.VariableNames;

%% Reading some other parameters as an example

% choose some specific time point
timeofinterest = 4; % seconds
time = table2array(data(:,"TestTime"));
TOI = find(abs(time - timeofinterest) <= 0.1); % finds point of interest based on time of interest 
TOI = TOI(1);

 % temperature from thermocouple(s)
temp = data(TOI,contains(name_list,'MCPtemp'));
temp = table2array(temp);

 % PWM voltage signal to the SMA. This is our input. Between 0 and 1.
pwm_data = data(TOI,contains(name_list,'pwm_duty_limb'));
pwm_data = table2array(pwm_data);

%% Now for the limb fiducial tags

% plot at the start of the test
tag_position_start = data(1,contains(name_list,'tag'));

% convert from table to arrays 
tag_position_array_start =table2array(tag_position_start);
tag_array_start(1,:) = tag_position_array_start(1:2:end);
tag_array_start(2,:) = tag_position_array_start(2:2:end);


%% Convert pixel to meter

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

xf = PixPos(1,:);
yf = PixPos(2,:);

%% PLOTTING
picturewidth_singlecolumn = 20;
font_size = 20;
line_width = 3;

color_map = parula(8);
color_map = [color_map(1,:);color_map(3,:);color_map(5,:);color_map(7,:)];
hold on
hfig= figure(1);  % save the figure handle in a variable
plot(xf, yf,'d-k',MarkerSize=10);

xlabel('$x\ (m)$');
ylabel('$y\ (m)$');
axis equal
axis([0 0.2 -0.12 0.06])

set(findall(hfig,'-property','FontSize'),'FontSize',font_size) % adjust fontsize to your document
set(findall(hfig,'-property','LineWidth'),'LineWidth',line_width)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

set(hfig,'Units','centimeters','Position',[3 3 picturewidth_singlecolumn picturewidth_singlecolumn]);


function PixPos = PixelToPosition(tag_array,convert)
    PixPos = tag_array*convert;
end 