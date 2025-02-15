# proj_td500_tempDynamics
MATLAB code for Temporal Dynamics of Object Recognition (RSVP500) project.
Datasets from Bell et al., 2011 and Morin et al., 2014
Based on RSVP tasks.

**rsvp500_massiveMatrix.data structure**
01 Neuron Number
02 trial_ID (1) - Condition Number (see below)
03 trial_ID (2) - Category (see below)
04 trial_m_baseline
05 trial_m_epoch1
06 spden_trial
07 MONKEY (need to add)
08-20) empty
21 spden_trial (-500 to +750 ms surrounding stimulus onset)

**Condition Numbers:**
hmiconfig.faces500  = 1:20; % condition numbers corresponding to face stimuli
hmiconfig.fruit500  = 21:40;
hmiconfig.places500 = 41:60;
hmiconfig.bodyp500  = 61:80;
hmiconfig.objct500  = 81:100;


% faces570_massiveMatrix.data structure
% 01) Neuron Number
% 02) condition number
% 03) expression?
% 04) monkeyface ID? 
% 05) gaze?
% 06) monkey?
% 07) 1, 2, or 3???
% 08) trial number?
% 09)
% 10)
% 11)


**faces570_massiveMatrix.data structure**
01 Neuron Number
02 trial_ID (1) - Stimulus Number (Condition Number, 1-88; see below)
03 trial_ID (2) - Expression (1=neutral,2=threat,3=fear,0=object/bodypart)
04 trial_ID (3) - Identity (1-8 identity 1-8, 0=object/bodypart)
05 trial_ID (4) - Gaze Direction (1=directed,2=averted,0=object/bodypart)
06 trial_ID (5) - Category (1=face,2=object,3=bodypart)
07 trial_ID (6) - Stimulus repetition
08 trial_ID (7) - Category repetition (faces, objects, body-parts)
09 trial_ID (8) - Identity repetition (face 1-8)
10 trial_ID (9) - Expression repetition (1-3)
11 trial_ID (10) - Gaze Direction repetition (1 or 2)
12 trial_ID (11) - Stimulus repetition (Correct Only)
13 trial_ID (12) - Category repetition (faces, objects, body-parts)
14 trial_ID (13) - Identity repetition (face 1-8)
15 trial_ID (14) - Expression repetition (1-3)
16 trial_ID (15) - Gaze Direction repetition (1 or 2)
17 NaN
18 ???
19 trial_m_epoch1
20 spden_trial (5000+-500:750)




% 21) spden_trial (5000+-500:750)

%Your monkey 1: female, mid-age
%Your monkey 2: female, mid-age
%Your monkey 3: male, mid-age
%Your monkey 4: male, mid-age
%Your monkey 5: female, young
%Your monkey 6: female, mid-ge
%Your monkey 7: female, can't tell the age
%Your monkey 8: female, 20 years old
% faces570
hmiconfig.facesND570 = [1 7 13 19 25 31 37 43];
hmiconfig.facesNA570 = [2 8 14 20 26 32 38 44];
hmiconfig.facesTD570 = [3 9 15 21 27 33 39 45];
hmiconfig.facesTA570 = [4 10 16 22 28 34 40 46];
hmiconfig.facesFD570 = [5 11 17 23 29 35 41 47];
hmiconfig.facesFA570 = [6 12 18 24 30 36 42 48];
hmiconfig.bodyp570 = 49:68;
hmiconfig.objct570 = 69:88;

hmiconfig.f570_neutral = [1 7 13 19 25 31 37 43 2 8 14 20 26 32 38 44];
hmiconfig.f570_threat = [3 9 15 21 27 33 39 45 4 10 16 22 28 34 40 46];
hmiconfig.f570_feargrin = [5 11 17 23 29 35 41 47 6 12 18 24 30 36 42 48];
hmiconfig.f570_directed =[1 7 13 19 25 31 37 43 3 9 15 21 27 33 39 45 5 11 17 23 29 35 41 47];
hmiconfig.f570_averted = [2 8 14 20 26 32 38 44 4 10 16 22 28 34 40 46 6 12 18 24 30 36 42 48];




% First Attempt
% Set Working directory
currentDirectory = pwd;
disp(['Current working directory:', currentDirectory]);
% Graph
dat = faces570_massiveMatrix.data;
plot(mean(dat(dat(:,3)==1, 21:end)))
plot(mean(dat(dat(:,1)==1, 21:end)))
plot(mean(dat(dat(:,1)==1 & dat(:,3)==1, 21:end)))
plot(-500:750,mean(dat(dat(:,1)==1 & dat(:,3)==1, 21:end)))
xlabel('time (ms)')
ylabel('sp/s')
hold on; % this allows you to paste multiple lines on the same figure
plot(-500:750,mean(dat(dat(:,1)==1 & dat(:,3)==2, 21:end)),'b-')
plot(-500:750,mean(dat(dat(:,1)==1 & dat(:,3)==2, 21:end)),'r-')
plot(-500:750,mean(dat(dat(:,1)==1 & dat(:,3)==3, 21:end)),'g-')
% Load data (assuming 'rsvp500_massiveMatrix.data' is already loaded in your workspace)
dat = faces570_massiveMatrix.data;
% Define the category and neuron of interest
neuronNumber = 3; % Change this to your desired neuron number
categoryIndex = 1; % Change this to the category you want to analyze (e.g., 1 for faces, 2 for fruits, etc.)
gazeConditions = unique(dat(:, 3)); % Get unique gaze conditions for the specified category
% Prepare the figure
figure;
hold on;
% Loop through each gaze condition and plot the mean responses
for gaze = gazeConditions'
    % Get the mean responses for the specified neuron and gaze condition
    meanResponse = mean(dat(dat(:, 1) == neuronNumber & dat(:, 3) == gaze, 21:end));
    % Time vector for x-axis
    timeVector = -500:750;
    % Plot the mean response
    plot(timeVector, meanResponse, 'DisplayName', ['Gaze ' num2str(gaze)]);
end
% Labels and legend
xlabel('Time (ms)');
ylabel('sp/s');
title(['Neuron ' num2str(neuronNumber) ' Responses for Category ' num2str(categoryIndex)]);
legend('show'); % Show legend for the gaze conditions
grid on;
hold off; % Release the plot hold
% Save the figure
saveas(gcf, ['neuron_' num2str(neuronNumber) '_category_' num2str(categoryIndex) '_responses.png']);
 