%% Code here: The path on line 19 specific to this subject folder change as needed
%
%directory of the file
clc
clear
close all

%% Code below: pick which ever one is randomized
%DJ_N_  6+7
%DJ_I_
%DJ_O_

%RDJ_N_
%RDJ_I_
%RDJ_O_

%Get the files based on the directory that we fed to function uigetfile 

[file_list, path_n] = uigetfile('H:\Thesis subject data\Sub_2\Session I\DJ_N_*.c3d', 'Select files','MultiSelect', 'on');

if iscell(file_list) == 0
    file_list = {file_list};
end

%Get the number of files selected - VERY IMPORTANT
file_list_number=length(file_list);

%% Pull out the data from the c3d files using btk 

for i=1:length(file_list)
    filename(i)=btkReadAcquisition(char(file_list(:,i)));

    [analog_data(i),analogsInfo(i)]= btkGetAnalogs(filename(i));
end

%% FORCES FIRST

% Filter the forces and then combine forces from force plate 1 and 2

fc =15; %cutoff frequency [Hz]
fs=1000; %sample frequency [Hz]
n_order = 2; % filter order

for i=1:length(analog_data)
    
[b,a]=butter(n_order, fc/(fs/2));

Left_force_z_filt=filtfilt(b,a,analog_data(i).Force_Fz1);
analog_data(i).Force_Fz1=Left_force_z_filt;

Right_force_z_filt=filtfilt(b,a,analog_data(i).Force_Fz2);
analog_data(i).Force_Fz2=Right_force_z_filt;

end

%% Combine force from FP1 and FP2 since data set has two force platforms L (Fz1) and R (Fz2) foot on separate force plate

%Also make cells for the time of each jump

DJ_forces{file_list_number,1}.jump=[];  %FZ for the entire cycle
time_raw{file_list_number, 1}.raw=[]; %time for the entire cycle in seconds
dj_jump=[];

for i = 1 :file_list_number
    
    dj_jump= analog_data(i).Force_Fz1 + analog_data(i).Force_Fz2;
    t_raw=0:1/fs:(length(dj_jump)-1)/fs;
    
    DJ_forces{i,1}.jump=-(dj_jump);
    time_raw{i, 1}.raw=t_raw';
    
end

%% Clear variables we no longer need to keep workspace clean

clearvars path_n Right_force_z_filt Left_force_z_filt n_order filename dj_jump b a  

% 1. Calculate participants BW using the final one second of one trial. This is for SUB_2 
%COM Velocity - Divide Vertical force (minus BW) by body mass and then
%integrating the product using trapezoid rule

% Access the jump forces for the desired participant

for i=1:file_list_number
    jump_velo = DJ_forces{i, 1}.jump;
    
    % Extract the force data for the final one second
    sampling_rate = fs;  % Sampling rate in samples per second
    final_second_force = jump_velo(end - sampling_rate + 1:end);
    
    % Calculate the average force during the final one second
    BW(i) = mean(final_second_force);
end

%% 2. Gravity, box height, and subject's mass

g=-9.81;
bh=0.3;

mass=abs(BW/g);



%% Get the foot strike and toe off events and save them in a structure for the left and right foot
% 
foot_strk=[];
toe_off=[];
%Left foot strike and toe off frames
for k=1:file_list_number
    for i = 1:length(DJ_forces{k, 1}.jump)
        
        %break out of the loop after the last element
        if i == length(DJ_forces{k, 1}.jump)
            break
        end
        %when it goes from less than 15N to more than 15N we have a foot
        %strike
        if DJ_forces{k, 1}.jump(i) < 15 && DJ_forces{k, 1}.jump(i+1)>15
            foot_strk(k,end+1)=i+1;
        end
        %when it goes from greater than 15N to less than 15N we have a toe
        %off
        if DJ_forces{k, 1}.jump(i) >15 && DJ_forces{k, 1}.jump(i+1)<15
           toe_off(k,end+1)=i;
        end
        
    end
end

%% Push all the zero values out to one side then keep only the first foot strike

[m,n]=size(foot_strk); 

for i=1:m
  tmp=nonzeros(foot_strk(i,:)); 
  foot_strk(i,:)=[tmp.', zeros(1,n-numel(tmp))]; 
end


%Same for the toe offs left
[m,n]=size(toe_off); 

for i=1:m
  tmp=nonzeros(toe_off(i,:)); 
  toe_off(i,:)=[tmp.', zeros(1,n-numel(tmp))]; 
end

%%
%Keep the first foot strike and toe off - discard everything else

foot_strk= foot_strk(:,1:2)-1;
toe_off=toe_off(:,1:2);


%% Frames for each jump and parsed time for 1. from contact to end of trial and, 2. support/contact phase


time{file_list_number, 1}.jump=[];

for i=1:file_list_number
    temp=time_raw{i, 1}.raw(foot_strk(i,1):toe_off(i,1));
    temp=temp-temp(1);
    time_contact{i, 1}.raw=temp';
    
    
    temp=time_raw{i, 1}.raw(foot_strk(i,1):end);
    temp=temp-temp(1);
    time{i, 1}.jump=temp';
    
end

% time{file_list_number, 1}.FZ=[]; this is the time from foot strike to the
% end of the trial 


% time_contact{i, 1}.raw is the time from foot strike to the toe off


% Clear variables we no longer need to keep workspace clean

clearvars tmp temp sup n k m final_second_force


%% Create a support phase (GC to TO) with the GRF Forces (FZ) and velocity
% 
support{file_list_number, 1}.FZ=[];
support_inseconds=[];
sup=[];



for i=1:file_list_number
    support_inseconds(i,1)=(toe_off(i,1) - foot_strk(i,1))/fs; %1000HZ is sampling frequency - output is in seconds - jump time
    second_landing(i,1)=(foot_strk(i,2) - foot_strk(i,1))/fs; 
    support{i, 1}.FZ=DJ_forces{i, 1}.jump(foot_strk(i,1):toe_off(i,1)); %save here the FZ frames between foot strike and toe off incluind the foot strike frame
    
end


%% 3. Acceleration

% a_forces{file_list_number,1}.jump=[]; %acceleration of the body during the entire jump
% adj_jump=[];

%Acceleration for the entire jump
for i = 1 :file_list_number
    
    adj_jump= DJ_forces{i,1}.jump(foot_strk(i):end);
    
    adj_jump=(adj_jump - BW(i))/mass(i);
    a_forces{i,1}.jump=adj_jump;
    adj_jump=[];
end



%COM velocity determined by dividing vertical force (minus BW) by body mass

%%

% %Acceleration for the support phase jump
% accel_sup{file_list_number,1}.jump=[];

for i = 1 :file_list_number
    
    adj_jump= support{i, 1}.FZ;
    
    adj_jump=(adj_jump -  BW(i))/mass(i);
    accel_sup{i,1}.jump=adj_jump;
    adj_jump=[];
end


clearvars adj_jump
%% Touch Down Velocity and height of COM standing

si=0; %height of COM standing - if you go below the standing height
%it goes negative position and if you go higher than 
%standing height it becomes positive


%Touchdown velocity was then estimated
% from box height based on the conservation of mechanical energy principle as the square
% root of 2 × 9.81 × box height (in m) 
vi=abs(2*g*bh).^0.5 *-1;

%% Velocity and Displacement

%Velocity for the entire trial starting with the footstrike 
v{file_list_number, 1}.jumps=[];

for j=1:file_list_number
    %for the support phase
    %velocity
    v{j, 1}.jumps=integ(time_raw{j, 1}.raw(foot_strk(j,1):end),a_forces{j,1}.jump,vi);

end

%% 1. Recalculate participants Velocity for the last second of one trial. This is for SUB_2
% 
% For exmaple if the mean velocity during the final one second of data collection 
% equalled −0.12 m·s−1 rather than 0 m·s−1 then this value was deducted from the touchdown 
% velocity that was estimated using the conservation of mechanical energy principle. Numerical integration of the force–time record obtained from the second FP was then performed again using this updated touchdown velocity value to generate “corrected” COM 
% velocity and displacement values throughout the entire data recording. DH was estimated 
% from the updated touchdown velocity as: touchdown velocity squared
% divided by 19.62 (i.e., 2 x gravitational acceleration).

end_v=[];

for i = 1: file_list_number
    jump_velo = v{i, 1}.jumps;

    % Extract the velocity data for the final one second
    final_second_velocity = jump_velo(end - sampling_rate + 1:end);
    % Calculate the average velocity during the final one second
    end_v(i,1)=mean(final_second_velocity);
    
end


%Now we want to subtract the mean velocity values from vi values and save
%this to a new variable upd_vi so we can recalculate the actual velocity
%and position
%5. Subtract the value of the last minute BW measurement for the velocity curve from the estimated touch down

Upd_vi=vi-end_v; %Value becomes more positive when we subtract it by negative numbers

%% 3. velocity then was estimated using the conservation of mechanical energy principle. Numerical integration
% Recalculate Velocity and Displacement

%Velocity for the support phase
v{file_list_number, 1}.jumps=[];
v_sup{file_list_number, 1}.jump=[];
%Displacement for the support phase
s{file_list_number, 1}.jumps=[];


for j=1:file_list_number
    %for the support phase
    %velocity
    v{j, 1}.jumps=integ(time_raw{j, 1}.raw(foot_strk(j,1):end),a_forces{j,1}.jump,Upd_vi(j));
    v_sup{j, 1}.jump=integ(time_raw{j, 1}.raw(foot_strk(j,1):toe_off(j,1)),accel_sup{j,1}.jump,Upd_vi(j));

    %position
    s{j, 1}.jumps=integ(time_raw{j, 1}.raw(foot_strk(j,1):end),v{j, 1}.jumps,si);
    s_sup{j, 1}.jump=integ(time_raw{j, 1}.raw(foot_strk(j,1):toe_off(j,1)),v_sup{j, 1}.jump,si);

end



%% Find where vertical velocity is 0 to differentiate between the braking and propulsion phases

%at what frame does velocity go from negative to positive

% Find the index where the velocity transitions from negative to positive

for i = 1:file_list_number
    
    velocity_data = v{i, 1}.jumps;

    %transitioning index is the frame/column for the entire velocity curve
    %starting with the footstrike to the end of the trial
    transition_index(i,1) = find(diff(sign(velocity_data)) > 0, 1); %the last negative value before it changes up to positive

    time_tr_ind(i,1)=time_contact{i, 1}.raw(transition_index(i,1));

end
clearvars velocity_data

%%

for i = 1: file_list_number
%     transition_index(i,1) = find(v_sup{i, 1}.jump >0,1,'first'); % +1 because first vlaue bigger than 0 
     trans_vsup_index(i,1)=transition_index(i,1)+1; %  v + 1 because first value larger than 0 
    
    
     temp= time_contact{i, 1}.raw(1:trans_vsup_index(i,1)-1); % needed for area function in figure
     x_help_ecc{i,1}.ecc=temp;
     
     temp = time_contact{i, 1}.raw(trans_vsup_index(i,1)-1:end);
     x_help_con{i,1}.con=temp;
 % % %     x_help_con = start_con:1:take_off_2fp; % needed for area function in figure
end

%% Plot data
figure(1);
% tiledlayout(1,2, 'Padding','tight','TileSpacing','tight')
t = tiledlayout(2,2,'Padding','tight','TileSpacing','tight');
S=gcf;
set(gcf,'WindowState','maximized','Color', [1 1 1]);
%Plot 1
nexttile
subtitle("(a) Vertical displacement and vertical velocity vs time",'FontSize',12)
yyaxis left

plot(time{2, 1}.jump,s{2, 1}.jumps, "b-",LineWidth=2)
hold on
ylabel('Position (m)')
xlabel('Time (s)')
ylim([-0.5 0.5])
yticks([-0.5 -0.25 0 0.25 0.5])

yyaxis right
plot(time{2, 1}.jump,v{2, 1}.jumps,"r-",LineWidth=2)
ylabel('Velocity (m/s)')
yline(0,'k-','Fontsize',7)

xlim([0 2])
xticks([0 0.5 1 1.5 2])

%labels
% lab=["Position", "Velocity"];
% legend(lab)
xline(time_tr_ind(2),'k--','Fontsize',7)
xline(support_inseconds(2),'g-.','Fontsize',7,LineWidth=1.3)
xline(second_landing(2),'r-.','Fontsize',7,LineWidth=1.3)


lab=["Position", "Velocity"];
legend(lab)
text(1.63,0.3,"Velocity during ",'Fontsize',6)
text(1.63,0.1,"the weighing period ",'Fontsize',6)
text(0.35,-2,"Flight time ",'Fontsize',12)

nexttile;

plot(time_contact{2, 1}.raw,v_sup{2, 1}.jump,"r-",LineWidth=2)

subtitle("(b) Velocity vs Time during the support phase",'FontSize',12)
ylabel('Velocity (m/s)')
xlabel('Time (s)')
yline(0,'k-','Fontsize',7)
xline(time_tr_ind(2),'k--','Fontsize',7)


nexttile

plot(time{2, 1}.jump,DJ_forces{2, 1}.jump(foot_strk(2):end), "k-",LineWidth=2)
subtitle("(c) Force vs Time",'FontSize',12)
ylabel('Ground Reaction Force (N)')
xlabel('Time (s)')
xlim([0 2])
xticks([0 0.5 1 1.5 2])
xline(toe_off(2),'g--','Fontsize',7)
text(1.63,679,"Force during ",'Fontsize',6)
text(1.63,590,"the weighing period ",'Fontsize',6)
xline(time_tr_ind(2),'k--','Fontsize',7)
xline(support_inseconds(2),'g-.','Fontsize',7,LineWidth=1.3)
xline(second_landing(2),'r-.','Fontsize',7,LineWidth=1.3)
yline(0,'k-','Fontsize',7)
text(0.35,1590,"Flight time ",'Fontsize',12)
xlim([0 2])
xticks([0 0.5 1 1.5 2])

nexttile
plot(time_contact{2, 1}.raw,support{2, 1}.FZ,"k-",LineWidth=2)
hold on
area(x_help_ecc{2, 1}.ecc,support{2, 1}.FZ(1:transition_index(2)),'FaceColor','r')
area(x_help_con{2, 1}.con,support{2, 1}.FZ(transition_index(2):end),'FaceColor','g')
subtitle("(d) Force vs Time during the support phase",'FontSize',12)
ylabel('Force (N)')
xlabel('Time (s)')
xline(time_tr_ind(2),'k--','Fontsize',7)
% area(x_help_con,Dj_sum_cut(start_con:take_off_2fp),'FaceColor','g')


%%
num=3;


figure(2);
% tiledlayout(1,2, 'Padding','tight','TileSpacing','tight')
t = tiledlayout(2,2,'Padding','tight','TileSpacing','tight');
S=gcf;
set(gcf,'WindowState','maximized','Color', [1 1 1]);
%Plot 1
nexttile
subtitle("(a) Vertical displacement and vertical velocity vs time",'FontSize',12)
yyaxis left

plot(time{num, 1}.jump,s{num, 1}.jumps, "b-",LineWidth=2)
hold on
ylabel('Position (m)')
xlabel('Time (s)')
ylim([-0.5 0.5])
yticks([-0.5 -0.25 0 0.25 0.5])

yyaxis right
plot(time{num, 1}.jump,v{num, 1}.jumps,"r-",LineWidth=2)
ylabel('Velocity (m/s)')
yline(0,'k-','Fontsize',7)

xlim([0 2])
xticks([0 0.5 1 1.5 2])

%labels
% lab=["Position", "Velocity"];
% legend(lab)
xline(time_tr_ind(num),'k--','Fontsize',7)
xline(support_inseconds(num),'g-.','Fontsize',7,LineWidth=1.3)
xline(second_landing(num),'r-.','Fontsize',7,LineWidth=1.3)


lab=["Position", "Velocity"];
legend(lab)
text(1.63,0.3,"Velocity during ",'Fontsize',6)
text(1.63,0.1,"the weighing period ",'Fontsize',6)
text(0.35,-2,"Flight time ",'Fontsize',12)

nexttile;

plot(time_contact{num, 1}.raw,v_sup{num, 1}.jump,"r-",LineWidth=2)

subtitle("(b) Velocity vs Time during the support phase",'FontSize',12)
ylabel('Velocity (m/s)')
xlabel('Time (s)')
yline(0,'k-','Fontsize',7)
xline(time_tr_ind(num),'k--','Fontsize',7)
xlim([0 time_contact{num, 1}.raw(end)])

nexttile

plot(time{num, 1}.jump,DJ_forces{num, 1}.jump(foot_strk(num):end), "k-",LineWidth=2)
subtitle("(c) Force vs Time",'FontSize',12)
ylabel('Ground Reaction Force (N)')
xlabel('Time (s)')
xlim([0 2])
xticks([0 0.5 1 1.5 2])
xline(toe_off(2),'g--','Fontsize',7)
text(1.63,679,"Force during ",'Fontsize',6)
text(1.63,590,"the weighing period ",'Fontsize',6)
xline(time_tr_ind(2),'k--','Fontsize',7)
xline(support_inseconds(num),'g-.','Fontsize',7,LineWidth=1.3)
xline(second_landing(num),'r-.','Fontsize',7,LineWidth=1.3)
yline(0,'k-','Fontsize',7)
text(0.35,1590,"Flight time ",'Fontsize',12)
xlim([0 2])
xticks([0 0.5 1 1.5 2])

nexttile
plot(time_contact{num, 1}.raw,support{num, 1}.FZ,"k-",LineWidth=2)
hold on
area(x_help_ecc{num, 1}.ecc,support{num, 1}.FZ(1:transition_index(num)),'FaceColor','r')
area(x_help_con{num, 1}.con,support{num, 1}.FZ(transition_index(num):end),'FaceColor','g')
subtitle("(d) Force vs Time during the support phase",'FontSize',12)
ylabel('Force (N)')
xlabel('Time (s)')
xline(time_tr_ind(num),'k--','Fontsize',7)
% area(x_help_con,Dj_sum_cut(start_con:take_off_2fp),'FaceColor','g')
xlim([0 time_contact{num, 1}.raw(end)])




%% CALCULATIONS OF VARIABLES

%% Recalculation of the drop height - corrected

drop_height=[];
for ii=1:file_list_number
    
    drop_height(ii)=((Upd_vi(ii).^2)/19.62);

end
%% Breaking phase time

brphase_time=time_tr_ind;

%% Propulsion phase time

pro_phase_time=[];

for i =1:file_list_number
    pro_phase_time(i,1)= x_help_con{i, 1}.con(end) - x_help_con{i, 1}.con(1);
end


%% Jump height - formula of touchdown. take off velocity squared divided by 19.62

jump_height=[];

for i =1:file_list_number
    temp= v_sup{i, 1}.jump(end)
    temp=((temp^2)/19.62);
    jump_height(i,1)=temp;
end

%%
for i =1:file_list_number
    xcon(i,1)=length(x_help_con{i, 1}.con);
    xecc(i,1)=length(x_help_ecc{i, 1}.ecc)-1;

end
%% Braking Displacement and Propulsion Displacement

break_displc=[];
prop_displc=[];

for i =1:file_list_number
    break_displc(i,1)= s_sup{i, 1}.jump(xecc(i,1)); %position starts at 0 so it is this frame/time/index 
    
    prop_displc(i,1)=-s_sup{i, 1}.jump(xecc(i,1))+ s_sup{i, 1}.jump(end);
end




%% Mean breaking and propulsion force


break_force=[];
prop_force=[];

for i =1:file_list_number

    temp=support{i, 1}.FZ(1:transition_index(i,1));
    break_force(i,1)=mean(temp);
    temp=support{i, 1}.FZ(transition_index(i,1)+1:end);
    prop_force(i,1)=mean(temp);

end

%% RSI

RSI=[];
mRSI=[];

for i =1:file_list_number

    flight_time(i,1)=second_landing(i,1)-support_inseconds(i,1);
    RSI(i,1)=flight_time(i,1)/support_inseconds(i,1);
    
    mRSI(i,1)=jump_height(i,1)/support_inseconds(i,1); % Velocity index
end



%% Turn the support ground reaction phase and displacement phase into a cycle
for i =1:file_list_number

    % Assuming you have your support phase force and displacement data stored in variables forceData and displacementData
    forceData= support{i, 1}.FZ;
    displacementData=s_sup{i, 1}.jump;
    % Define the original phase range (e.g., 0-100%)
    originalPhase = linspace(1, 101, length(forceData)); 

    % Define the desired phase range (e.g., 1-100%)
    desiredPhase = linspace(1, 101, 100);

    % Use spline interpolation to obtain the interpolated force and displacement values at the desired phase range
    Force_cycle(i,:) = spline(originalPhase, forceData, desiredPhase);
    displacement_cycle(i,:)= spline(originalPhase, displacementData, desiredPhase);
    
end

%% Impact displacement - percentage of peak displacement that was completed in the first 20% of ground contact time
Force_normed=[];
for i =1:file_list_number
    % forces scaled  to  body-weight  to  account  for  the  non-  linear  relationship  between  muscular strength and body size
    
    Force_normed(i,:)=Force_cycle(i,:)/BW(1,i);

end    

%%
for i =1:file_list_number
    
    %Determine the peak displacement during the support phase:
    
    peak_displacement(i,:) = min(displacement_cycle(i,:), [], 2); % Find the maximum value along each trial

    %Determine the duration of the first 20% of the ground contact time:
    ground_contact_time = 100; % Assuming the ground contact time is represented by 100 data points
    duration_20_percent = ground_contact_time * 0.2; % Calculate 20% of the ground contact time
    
    %Calculate the displacement covered during the first 20% of the ground contact time:
    
    displacement_20_percent(i,:) = displacement_cycle(i, 1:duration_20_percent); % Select the corresponding portion of the displacement data

    %Calculate the impact displacement as a percentage of the peak displacement:
    impact_displacement(i,:) = (displacement_20_percent(i,:) ./ peak_displacement(i,:)) * 100; % Calculate the percentage of peak displacement
end
     


%% Correlation to ensure the jump was spring like


correlation_matrix = corrcoef(Force_normed(2,:), displacement_cycle(2,:));
correlation_coefficient = correlation_matrix(1, 2);

disp(['Correlation coefficient: ', num2str(correlation_coefficient)]);



%% Variables of interest - add them to the a table 


% %%
for i =1:file_list_number
    
    Performance_table.DJ_O(1,i)=(drop_height(i));
    Performance_table.DJ_O(2,i)=-(Upd_vi(i));
    Performance_table.DJ_O(3,i)=(jump_height(i));
    Performance_table.DJ_O(4,i)=(brphase_time(i));
    Performance_table.DJ_O(5,i)=(pro_phase_time(i));
    Performance_table.DJ_O(6,i)=(break_displc(i));
    Performance_table.DJ_O(7,i)=(prop_displc(i));
    Performance_table.DJ_O(8,i)=(break_force(i));
    Performance_table.DJ_O(9,i)=(prop_force(i));
    Performance_table.DJ_O(10,i)=(RSI(i));
    Performance_table.DJ_O(11,i)=(mRSI(i));

end
