%Analysis
clc
clear all 
close all


%% Import data from spreadsheet


%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 29);

% Specify sheet and range
opts.Sheet = "MetasoftStudio";
opts.DataRange = "A3:AC983";

% Specify column names and types
opts.VariableNames = ["t", "Phase", "Marker", "WR", "VO2", "VO2kg", "VCO2", "VO2HR", "RER", "VEVO2", "VEVCO2", "VE", "VT", "BF", "HR", "SVestVO2max", "QestVO2max", "CI", "EDVest", "EDFR", "LCWi", "SVR", "SVRi", "EE", "EEFAT", "EECHO", "EEkg", "CHO", "FAT"];
opts.VariableTypes = ["char", "categorical", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "categorical", "categorical", "categorical", "categorical", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, ["t", "Marker"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["t", "Phase", "Marker", "EDVest", "EDFR", "LCWi", "SVR", "SVRi"], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable("D:\Stuff\CPET Rebecca\CPET.xlsx", opts, "UseExcel", false);

%% Convert to output type
t = tbl.t;
Phase = tbl.Phase;
Marker = tbl.Marker;
WR = tbl.WR;
VO2 = tbl.VO2;
VO2kg = tbl.VO2kg;
VCO2 = tbl.VCO2;
VO2HR = tbl.VO2HR;
RER = tbl.RER;
VEVO2 = tbl.VEVO2;
VEVCO2 = tbl.VEVCO2;
VE = tbl.VE;
VT = tbl.VT;
BF = tbl.BF;
HR = tbl.HR;
SVestVO2max = tbl.SVestVO2max;
QestVO2max = tbl.QestVO2max;
CI = tbl.CI;
EDVest = tbl.EDVest;
EDFR = tbl.EDFR;
LCWi = tbl.LCWi;
SVR = tbl.SVR;
SVRi = tbl.SVRi;
EE = tbl.EE;
EEFAT = tbl.EEFAT;
EECHO = tbl.EECHO;
EEkg = tbl.EEkg;
CHO = tbl.CHO;
FAT = tbl.FAT;

%%

t = split(t, ':');

m=zeros(size(t,1),size(t,2));
m=str2double(t);
m=m(:,2:3);

joined = cell(size(m, 1), 1);
for i = 1:size(m, 1)
    joined{i} = sprintf('%d.%d', m(i, 1), m(i, 2));
end
time = cellfun(@str2double, joined);
% %% Ventilatory response panels

[mvv, mvvid]=max(VE);
mvv2=time(mvvid,:);
BR = (mvv - VE) ./ mvv;
br=BR*100;


%%
% Information on the cardiovascular response and oxygen transport is
% reflected in panels 1 → 2 → 3 (original version: panels
% 3 → 2 → 5).
% Information on pulmonary gas exchange and V/Q mis-
% match can be found in panels 4 → 6 → 7 (original version:
% panels 6 → 4 → 9). A possible limitation of ventilatory
% capacity is shown in panels 5 → 8 → 9 (original version:
% panels 1 → 8 → 7).



%% MOV AVE
windowSize = 30;  % Window size in seconds
VE_average = movmean(VE, windowSize, 'Endpoints', 'discard');
t_average = movmean(time, windowSize, 'Endpoints', 'discard');

RER_average= movmean(RER, windowSize, 'Endpoints', 'discard');

BR_average= movmean(br, windowSize, 'Endpoints', 'discard');


VT_average= movmean(VT, windowSize, 'Endpoints', 'discard');

BF_average= movmean(BF, windowSize, 'Endpoints', 'discard');

CI_average=movmean(CI, windowSize, 'Endpoints', 'discard');
%%

threshold = 20;  % Minimum value threshold
duration = 10;  % Minimum duration in subsequent values

index = 0;
for i = 1:length(BF_average)-duration
    if all(BF_average(i:i+duration-1) > threshold)
        index = i;
        break;
    end
end

if index > 0
    disp(['Found at index: ' num2str(index)]);
else
    disp('No such value found.');
end
%% BF_20
 
x_point = [0, VE_average(1,1)];  % x-coordinates of the line endpoints
y_point = [0, VT_average(1,1)];  % y-coordinates of the line endpoints

slope = VT_average(1,1) / VE_average(1,1);%Determine the slope of the line. The slope can be calculated as:
x = [0, VE_average(1,1), VE_average(1,1) + 100];  % Adjust the range as needed
y = slope * x;%Define a range of x-values that extends beyond VE_average(1,1) (to represent the line extending into infinity). For example, you can use:
%% BF_50

% Initialize variables
threshold2 = 50;     % Threshold value
numConsecutive2 = 10;  % Number of consecutive values required

% Find the index of the first value that is over the threshold
startIndex = find(BF_average > threshold2, 1);

% Check if startIndex is valid and there are at least numConsecutive values after it
if ~isempty(startIndex) && (startIndex + numConsecutive2 - 1) <= length(BF_average)
    % Check if the consecutive values meet the condition
    if all(BF_average(startIndex : startIndex + numConsecutive2 - 1) > threshold2)
        % Found a value that is over 50 and stays over 50 for the next 10 values
        requiredValue = BF_average(startIndex);
        disp(['The value ', num2str(requiredValue), ' meets the condition.']);
    else
        disp('No value meets the condition.');
    end
else
    disp('No value meets the condition.');
end
%%
x_point2 = [0, VE_average(startIndex,1)];  % x-coordinates of the line endpoints
y_point2 = [0, VT_average(startIndex,1)];  % y-coordinates of the line endpoints

slope2 = VT_average(startIndex,1) / VE_average(startIndex,1);%Determine the slope of the line. The slope can be calculated as:
x2 = [0, VE_average(startIndex,1), VE_average(startIndex,1) + 100];  % Adjust the range as needed
y2 = slope2 * x2;%Define a range of x-values that extends beyond VE_average(1,1) (to represent the line extending into infinity). For example, you can use:
%% % Shade the region between the two BF lines 
% Shade the region between the lines
x_fill = [x, fliplr(x2)];  % Combine x-coordinates in the correct order for filling
y_fill = [y, fliplr(y2)];  % Combine y-coordinates in the correct order for filling

%%  % Maximal voluntary ventilation value
threshold3 = 0.85 * mvv;  % Threshold value (85% of mvv)

% Find the indices where VE/mvv is greater than the threshold
indices = find((VE_average) > threshold3);
indices2=VE_average(indices(1,1),1);
% Display the indices
disp(indices);


%%
figure(1);
tiledlayout(1,3, 'Padding','tight','TileSpacing','tight')

%Plot 1

nexttile
plot(time,VE,"r.")
plot(t_average,VE_average,"r.")
a = get(gca,'XTickLabel');
set(gca,'FontName','Times','Fontsize',7)

subtitle("Minute ventilation (̇VE) and work rate (WR) vs. time")
xl=xline(mvv2,'k--',{'E',});
xl.LabelHorizontalAlignment = 'center';
xl.LabelVerticalAlignment = 'middle';

yyaxis left
ylim([0 200])
xl=xline(3,'k--',{'B',});
xl.LabelHorizontalAlignment = 'center';
xl.LabelVerticalAlignment = 'middle';
yline(mvv,'k--',{'MVV'},'FontSize',4)
ylabel("VE [L/min]",'Color','r')
yyaxis right
plot(time,WR,"ks")
xlabel("Time [minutes]")
ylabel("WR [W]",'Color','k')
yticks([50 80 110 140 170 200])
h=gca;
h.FontWeight='light';

%Plot 2
nexttile
plot(t_average,RER_average,"r.")
title("Ventilatory response panels")
subtitle("Respiratory exchange rate (RER) and breathing reserve (BR)")
a = get(gca,'XTickLabel');
set(gca,'FontName','Times','Fontsize',7)
xl=xline(mvv2,'k--',{'E',});
xl.LabelHorizontalAlignment = 'center';
xl.LabelVerticalAlignment = 'middle';
ylim([0.75 1.15])
yline(0.8,'k--',{'BR'},'fontsize',7)
yline(1.0,'k--',{'RER = 1'},'fontsize',7)
yyaxis left
yline(mvv,'k--',{'E'})
ylabel("RER ",'FontWeight','bold','Color','r','fontsize',7)
xl=xline(3,'k--',{'B',});
xl.LabelHorizontalAlignment = 'center';
xl.LabelVerticalAlignment = 'middle';
yyaxis right
plot(t_average,BR_average,"b.")
xlabel("Time [minutes]")
ylabel("BR [%]",'Color','b')
h=gca;
h.FontWeight='light';


%Plot 3 VT v CE
nexttile
plot(VE_average,VT_average,"r.")
a = get(gca,'XTickLabel');
set(gca,'FontName','Times','Fontsize',7)
hold on
ylim([0 4])
ylabel("VT [L]",'FontSize',7,'Color','r')

subtitle("Tidal volume (VT), minute ventilation (̇VE) and breathing frequency (BF)",'FontSize',7)
xlabel("̇VE [L/minutes]")
% plot(x, y, "--","Color", [0.5 0.5 0.5]);  % Adjust the line color and style as desired  % Extend the line through the graph
% plot(x2, y2, "--","Color", [0.5 0.5 0.5]);  % Adjust the line color and style as desired  % Extend the line through the graph
xline(indices2,'k-',{'V E/MVV >  85%',},'FontSize',7);
xlim([0 150])
h=fill(x_fill, y_fill, [0.6 0.6 0.6],'FaceAlpha', 0.2);  % Adjust the fill color and transparency as desired
% plot(x, y, "--","Color", [0.5 0.5 0.5]);  % Adjust the line color and style as desired  % Extend the line through the graph
% plot(x2, y2, "--","Color", [0.5 0.5 0.5]);  % Adjust the line color and style as desired  % Extend the line through the graph
h.EdgeColor="none";

text(43,2.5,'BF=20','fontsize',7)
text(135,2.7,'BF=50','fontsize',7)

S=gcf;
set(gcf,"WindowState","maximized","Color",[1 1 1])

h=gca;
h.FontWeight='light';



%% Cardiovascular panels


vo2peak=2.92;
VCO2_average = movmean(VCO2, windowSize, 'Endpoints', 'discard');
VO2_average = movmean(VO2, windowSize, 'Endpoints', 'discard');
WR_average= movmean(WR, windowSize, 'Endpoints', 'discard');
% Calculate the change in VO2 and WR
Aerobic_contribution= (VO2_average./WR_average);
Aerobic_contribution=Aerobic_contribution*1000;

[peakVO2,indpkvo2]=max(VO2_average);
WR_peak=WR_average(indpkvo2,1);

aerob_max=peakVO2/WR_peak;
aerob_max=aerob_max*1000;
VO2kg_average = movmean(VO2kg, windowSize, 'Endpoints', 'discard');
new_windSize= windowSize*6;
tim3aver=movmean(t_average, new_windSize, 'Endpoints', 'discard');

HR_average = movmean(HR, windowSize, 'Endpoints', 'discard');

VO2HR_average=movmean(VO2HR, windowSize, 'Endpoints', 'discard');

VO2hrpk=VO2HR_average(indpkvo2,1);

%%

aero_3minaverage = movmean(Aerobic_contribution, new_windSize, 'Endpoints', 'discard');
new_windSize3=60;
avg_aero=mean(aero_3minaverage);
VO2_180average = movmean(VO2, new_windSize, 'Endpoints', 'discard');
HR_180average = movmean(HR, new_windSize, 'Endpoints', 'discard');

VO2_60 = movmean(VO2, 60, 'Endpoints', 'discard');
HR_60 = movmean(HR, 60, 'Endpoints', 'discard');

VT1=2.02;%%


%% Cardiovascular panels
figure(2);
tiledlayout(1,3, 'Padding','tight','TileSpacing','tight')

%Plot 1
nexttile

plot(t_average,VO2_average,"b.")
subtitle("O2 uptake and CO2 output  vs. time and peak VO2 vs. work rate (WR)",'Fontsize',7)
a = get(gca,'XTickLabel');
set(gca,'FontName','Times','Fontsize',7)
xl=xline(mvv2,'k--',{'E',},'Fontsize',7);
xl.LabelHorizontalAlignment = 'center';
xl.LabelVerticalAlignment = 'middle';
xl=xline(3,'k--',{'B',},'Fontsize',7);
xl.LabelHorizontalAlignment = 'center';
xl.LabelVerticalAlignment = 'middle';

yyaxis left
ylim([0 4])
yline(vo2peak,"k--")
xlabel("Time [minutes]")
ylabel("VO˙2 [L/min]",'Fontsize',7,'Color','b')

yyaxis right

plot(t_average,VCO2_average,"r.")
text(10,0.3,'Average ΔV̇O2/ΔWR=13.3mL/min/W','Fontsize',7)
text(3,3,'Aerobic metabolism at peak VO2  ΔV̇O2/ΔWR=14.58mL/min/W','Fontsize',7)
ylabel("VCO˙2 [L/min]",'Fontsize',7,'Color','r')
ylim([0 4])
%%
%Plot 2
nexttile
title("Cardiovascular panels")
plot(t_average,HR_average,"r.")

text(14,114,'100% O2 pulse pred','Fontsize',7)
subtitle("Heart rate and oxygen pulse vs. time",'Fontsize',7)
a = get(gca,'XTickLabel');
set(gca,'FontName','Times','Fontsize',7)
yyaxis left
ylim([80 180])
ylabel("HR [beats/minute]",'Fontsize',7,'Color','r')
xlabel("Time [minutes]")
yyaxis right

plot(t_average,VO2HR_average,"b.")
ylabel("O2 pulse (V̇ O2/HR) [mL/beat]",'Fontsize',7,'Color','b')
ylim([0 50])

xl=xline(mvv2,'k--',{'E',},'Fontsize',7);
xl.LabelHorizontalAlignment = 'center';
xl.LabelVerticalAlignment = 'middle';
xl=xline(3,'k--',{'B',},'Fontsize',7);
xl.LabelHorizontalAlignment = 'center';
xl.LabelVerticalAlignment = 'middle';

yline(VO2hrpk,'k--','Fontsize',7)

title("Cardiovascular panels")

%Plot 3

nexttile
plot(VO2_average,VCO2_average,"r.")
xlabel("V̇O2[L/min]")
subtitle("O2 uptake vs CO2 output and HR vs VO˙2.",'Fontsize',7)
a = get(gca,'XTickLabel');
set(gca,'FontName','Times','Fontsize',7)
yyaxis left

xline(VT1,'k-',"LineWidth",2)
ylim([0 4])
ylabel("V̇ CO2 [L/min]",'Fontsize',7,'Color','r')
yyaxis right
plot(VO2_60,HR_60,"b.")
hold on
plot(VO2,HR,"b.")
ylim([0 250])
text(1.9,200,'AT','Fontsize',7)
ylabel("HR [beats /min]",'Fontsize',7,'Color','b')
% yticks([50 80 110 140 170 200])
% h=gca;
% h.FontWeight='light';



%%
vevo2=VE_average./VO2_average;

vevco2=VE_average./VCO2_average;

%% Find the time when t was hit the value of AT or VT1 and was greater for the next 8 frames


% Initialize variables
threshold2 = VT1;     % Threshold value
numConsecutive2 = 4;  % Number of consecutive values required

% Find the index of the first value that is over the threshold
startIndex = find(VO2_average > threshold2, 1);

% Check if startIndex is valid and there are at least numConsecutive values after it
if ~isempty(startIndex) && (startIndex + numConsecutive2 - 1) <= length(VO2_average)
    % Check if the consecutive values meet the condition
    if all(VO2_average(startIndex : startIndex + numConsecutive2 - 1) > threshold2)
        % Found a value that is over 50 and stays over 50 for the next 10 values
        requiredValue = VO2_average(startIndex);
        disp(['The value ', num2str(requiredValue), ' meets the condition.']);
    else
        disp('No value meets the condition.');
    end
else
    disp('No value meets the condition.');
end

%Now let's find where this is in in the time vector

t_vt1=t_average(startIndex,1);

%% Find the time when VCO2 hit the value of AT or VT1 and was greater for the next 8 frames

VCO2_vt1=1.94;

% Initialize variables
threshold2 = VCO2_vt1;     % Threshold value
numConsecutive2 = 2;  % Number of consecutive values required

% Find the index of the first value that is over the threshold
startIndex = find(VCO2_average > threshold2, 1);

% Check if startIndex is valid and there are at least numConsecutive values after it
if ~isempty(startIndex) && (startIndex + numConsecutive2 - 1) <= length(VCO2_average)
    % Check if the consecutive values meet the condition
    if all(VCO2_average(startIndex : startIndex + numConsecutive2 - 1) > threshold2)
        % Found a value that is over 50 and stays over 50 for the next 10 values
        VErequiredValue = VCO2_average(startIndex);
        disp(['The value ', num2str(requiredValue), ' meets the condition.']);
    else
        disp('No value meets the condition.');
    end
else
    disp('No value meets the condition.');
end


VE_vt1=VE_average(startIndex,1);

%% Squeezing in a regression line

%parse VE/VCO2 data and create a variable wherever the AT/VT1 line cuts
%thorugh the second plot

VE_vslop=VE_average(1:startIndex,1);

VCO2_vslop=VCO2_average(1:startIndex,1);




% Assuming x and y are your data points
degree = 1;  % Degree of the polynomial fit (1 for linear regression)
coefficients = polyfit(VE_vslop, VCO2_vslop, degree);
%Generate the x-values for the regression line using the range of your data
x_regression = linspace(min(VE_vslop), max(VE_vslop), 100);  % Adjust the number of points as needed
x_regression=x_regression';

%Evaluate the regression line using the polyval function.
y_regression = polyval(coefficients, x_regression);
y_regression=y_regression';



%%

figure(3);
tiledlayout(1,2, 'Padding','tight','TileSpacing','tight')


%Plot 1
nexttile

plot(t_average,vevo2,"b.","Linewidth",2)
title("Pulmonary gas exchange panels");
a = get(gca,'XTickLabel');
set(gca,'FontName','Times','Fontsize',7)
subtitle("Minute ventilation vs O2 uptake vs CO2 output as a function of time",'FontSize',7)
ylim([0 60])
xlabel("Time [minutes]")
yyaxis left
ylabel("EqO2 (~ V̇ E / V̇ O2)",'FontSize',7,'Color','b')
xl=xline(mvv2,'k--',{'E',},'FontSize',7);
xl.LabelHorizontalAlignment = 'center';
xl.LabelVerticalAlignment = 'middle';
xl=xline(3,'k--',{'B',},'FontSize',7);
xl.LabelHorizontalAlignment = 'center';
xl.LabelVerticalAlignment = 'middle';
xl=xline(t_vt1,'k-',{'AT/VT1',},'FontSize',7);

text(10,40,"Values at AT EqO2AT = 29; EqCO2 = 30.8  'normal'",'Fontsize',7)


yyaxis right
plot(t_average,vevco2,"r","Linewidth",2)
hold on
ylabel("EqCO2 (~ V̇ E / V̇ CO2)",'FontSize',7,'Color','r')
ylim([0 60])

%Plot 2
nexttile

plot(VCO2_average,VE_average,"r.","Linewidth",2)
hold on

xlabel("V̇CO2[L/min]")
ylabel("V̇E [L/min]",'FontSize',7,'Color','b')
a = get(gca,'XTickLabel');
xl=xline(VCO2_vt1,'k--',{'AT/VT1',},'FontSize',7);
set(gca,'FontName','Times','Fontsize',7)
subtitle("Ventilation (VE˙ ) and CO2 production(VCO˙2): VE˙/VCO˙2 slope",'FontSize',7)
ylim([0 150])
xlim([0 4])
yline(VE_vt1)

plot(y_regression,x_regression,"k-", "Linewidth",2);  % Regression line

text(0.5,62,"mean of VE˙/VCO˙2 slope)=> 30.7 'normal' ",'Fontsize',7)


%%
figure(4);

plot(t_average,Aerobic_contribution,'c',"LineWidth",3)
title("Contribution of aerobic metabolism to exercise (aerobic capacity)")

hold on
xlabel("Time [minutes]")
ylabel("∆V̇ O2 /∆WR [mL/min/W]",'FontSize',7,'Color','c',"LineWidth",3)

plot(tim3aver,aero_3minaverage,'r',"LineWidth",3)

% yticks([50 80 110 140 170 200])
% h=gca;
% h.FontWeight='light';

%%

