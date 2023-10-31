% Youtube video tutorial


clc
clear
close all

%% User enters trials for each one of the conditions


%Define the conditions and their order

conditions = {'DJ_N_','DJ_I_','DJ_O_','RDJ_N_','RDJ_I_','RDJ_O_'};

%Initialize the cell array for file names
fnames=cell(6,1);

for i =1:length(conditions)
    
    condition=conditions{i};
    
    %Prompt the user for the trial number for the current condition
    trial_number=input(['Enter the trial number for ' condition],'s');
    
    
    %Combine the condition and trial number to create the cell entry
    fnames{i}=[condition trial_number];
    
end


%%

%Define the base directory where the .c3d files are located
base_directory= 'C:\Users\trimm\Desktop\School related\YouT lectures\Video 3\Session I\';

%Initialize an empty cell array to store the selected file paths
selected_file_paths={};

for i = 1:length(fnames)

    file_pattern= [base_directory, fnames{i},'*.c3d'];
    [files, path_n]=uigetfile(file_pattern,['Select files for ', fnames{i}],'Multiselect','on');
    
    if iscell(files)

        selected_files_paths=[selected_file_paths, strcat(path_n, files)];
    else

        selected_file_paths{end+1}=fullfile(path_n, files);
    end
    %BTK Toolkit
    file_list{i}=files;

    filename(i)=btkReadAcquisition(char(file_list(:,i)));
    
    %parse data
    [angles(i), anglesInfo(i)]=btkGetAngles(filename(i));
%     [forces(i), forcesInfo(i)] = btkGetForces(filename(i));
    [moments(i).conds, momentsInfo(i)]=btkGetMoments(filename(i));
%     [moments(i).RKneeMoment__cgm23, momentsInfo]=btkGetMoments(filename(i));

%     [Events(i), eventsInfo(i)]=btkGetEvents(filename(i));
%     [forceplates{i}, forceplatesInfo{i}] = btkGetForcePlatforms(filename(i));
    [analog_data(i),analogsInfo(i)]= btkGetAnalogs(filename(i));
%     [ff(i)] = btkGetFirstFrame(filename(i));
    [markers(i).conds, markersInfo, markersResidual] = btkGetMarkers(filename(i));
        
end


file_list_number=length(file_list);

%%

subjects.male=struct();

subjects.female=struct();

while true
    
    path_n= input('Enter the path or "q" to quit: ','s');
    
    if strcmpi(path_n,'q')
        break;
    end
    
    
    path_parts = strsplit(path_n, '\');
    
    sub_index = find(contains(path_parts, 'Sub_'));
    
    Sub = path_parts{sub_index};
    
    sex=input('For female subjects, press (f). For male subjects, press(m):','s');
    
    subject_info.ID=Sub;
    
    if strcmpi(sex,'f')
        subjects.female.(Sub)=subject_info;
    elseif strcmpi(sex,'m')
        subjects.male.(Sub)=subject_info;
    else
        disp('Invalid input. Please enter "f" for female or "m" for male.');
    end
end









