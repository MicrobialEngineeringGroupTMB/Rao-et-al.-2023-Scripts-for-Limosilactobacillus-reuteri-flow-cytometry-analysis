%%%%% FCM_analysis_live_dead_damaged_fixed_gates %%%%%

% Written by Magnus Carlquist 2018-11-08
% Updated by Daniel P. Brink 2023-02-21 for the initial GitHub commit

% Purpose: this script was made to visualise FCM data and apply fixed gates
%          to distinguish three different subpopulations: intact cells, 
%          damaged cells A, and  damaged cells B. The script was developed 
%          using FCS3.1 formatted data from Limosilactobacillus reuteri 
%          stained with SYBRgreen and Propidium iodide and captured with a 
%          BD Accuri C6+ flow cytometer.  
%
% Calculations: number of cells, frequencies of intact, damagedB and 
%          damagedA cells.
%
% Plots:   number of cells/ml over time, light scatter bivarate plot 
%          (FSC-H, vs SSC-H) and viability bivariate plots (FL1-H vs FL3-H). 
%          Gated cells are visualised as green (intact), yellow (damagedA),
%          and red(damagedB).
%
% Dependencies: * FCSData_to_MData.m, which is a wrapper script for
%               fca_readfcs.m by Laszlo Balkay, available from from the Matlab Central
%               https://www.mathworks.com/matlabcentral/fileexchange/9608-fca_readfcs 
%               * Matlab toolboxes: Statistics and Machine Learning, 
%               Curve Fitting, Fixed-point Designer


%%%               %%%
%%% 1. User input %%%
%%%               %%%

%Specify the location of the fcs files. The folder is allowed to contain
%other file types, since the wrapper script is set to only parse *.fcs
%files.
source_dir='C:\Users\Daniel\Desktop\FCM_data';

%For the cell concentration estimates, the injection volume in µl and
%dilution factor of the samples needs to be specified by the user
%in the paper, 20 µl and 1000x dilution was used. Please change these
%values accordingly:
injection_volume=20;
dilution_factor=1000;

%%%                   %%%
%%% 2. Initialization %%%
%%%                   %%%

% The channel_labels variable contains the labels for the channels recorded
% in the FCS file. The BD Accuri C6+ files contains 14 channels. Example
% for the first loaded fcs file, FCMdata(1).
% Channel   Position in array 
% FSC-A     FCMdata(1).fcsdat(:,1)
% SSC-A     FCMdata(1).fcsdat(:,2)
% FL1-A     FCMdata(1).fcsdat(:,3)
% FL2-A     FCMdata(1).fcsdat(:,4)
% FL3-A     FCMdata(1).fcsdat(:,5)
% FL4-A     FCMdata(1).fcsdat(:,6)
% FSC-H     FCMdata(1).fcsdat(:,7)
% SSC-H     FCMdata(1).fcsdat(:,8)
% FL1-H     FCMdata(1).fcsdat(:,9)
% FL3-H     FCMdata(1).fcsdat(:,10)
% FL3-H     FCMdata(1).fcsdat(:,11)
% FL4-H     FCMdata(1).fcsdat(:,12)
% Width     FCMdata(1).fcsdat(:,13)
% Time      FCMdata(1).fcsdat(:,14)

% Key channels for this script are FSC-H (light scatter), 
% SSC-H (light scatter), FL1-H (SYBRgreen stain) and 
% FL3-H (Propidium iodide stain)

close all
clear timesincestart

%Load the fcs files
[FCMdata,number_of_fcs_files_in_folder] = FCSData_to_MData(source_dir);

%%%                   %%%
%%% 3. Calculations   %%%
%%%                   %%%

for z=1:number_of_fcs_files_in_folder 
    % Load the data
    % In the FCMdata file information of fcs-files include date and starttime
    % when they were recorded. This information can be used to fetch time
    % for sample points. "FCMdata(z).fcshrd.date", and "FCMdata(z).fcshrd.starttime".

    timesincestart(z) = datenum(FCMdata(z).fcshrd.date)-datenum(FCMdata(1).fcshrd.date)+datenum(FCMdata(z).fcshrd.starttime)-datenum(FCMdata(1).fcshrd.starttime);
    timesincestart(z) = timesincestart(z).*24; %convert the unit from days to hours 

    % Import the matrix for each sample into a temperary matrix called 'Data'
    clear Data; %to make sure the temporary matrix 'Data' is reinitated for each new sample
    Data=FCMdata(z).fcsdat(:,:,1);
           
    % Copy relevant channel number data from the Matrix 'Data' into vectors
    FSC_A=Data(:,1,:);
    SSC_A=Data(:,2,:);
    FSC_H=Data(:,7,:);
    SSC_H=Data(:,8,:);
    FL1_H=Data(:,9,:);
    FL3_H=Data(:,11,:);
    Width=Data(:,13,:);
           
    % Populate a matrix with log transformed data from the channels
    X = [log10(FSC_H), log10(SSC_H), log10(FL1_H), log10(FL3_H), log10(FSC_A), log10(SSC_A), Width];
    
    % Get the length of data matrix, before cleanup of the raw data
    L1(z)=length(X); 
    
    %% Clean up the data, based on a linear equation

    % Set new threshhold based on FL3-H X(:,4,:) and FL1-H X(:,3,:) to 
    % remove cell debris, unstained events, smaller events are not cells.
    
    % For the dataset used for developing the script, it was found that the
    % line y=-2.2202*x + 10.221 was suitable to remove such events from the
    % log transformed data
    % example code to visualize the line:
    % x_values=linspace(0,8)
    % y_values=-2.2202*x_values + 10.221
    % scatter(X(:,3),X(:,4), 3, 'k')
    % hold on
    % plot(x_values,y_values)

    % Remove all events from X that are smaller than or equal to the linear
    % Equation
    X(any(X(:,3,:)<=-2.2202*X(:,4,:) + 10.221,2),:) = []; % the 2 in any(M,2) makes the function test by row instead of by colum

    % Get the length of data matrix, after cleanup of the raw data
    L2(z)=length(X); 
    % Calculate the percentage of events that were removed during the cleanup. 
    percent_removed_rows(z)=(L1(z)-L2(z))/L1(z)*100; 
    
    % Calculate the mean values for the cleaned-up data for the key
    % channels
    mean_fsc_h(z)= mean(X(:,1,:));
    mean_ssc_h(z)= mean(X(:,2,:));
    mean_fsc_a(z)= mean(X(:,5,:));
    mean_ssc_a(z)= mean(X(:,6,:));

    %% Set fixed gates based on linear equations    
    % Set a fixed gate for the population with low FL1-H (SYBRgreen) signal.
    % This represents one the damaged cell subpopulations, here called
    % damaged B. For the dataset used for developing the script, it was 
    % found that the line y=0.908*x + 0.0265 was suitable to remove such 
    % events from the log transformed data.
    XdamagedB=X; 
    XdamagedB(any(XdamagedB(:,3,:)>=0.908*XdamagedB(:,4,:) + 0.0265,2),:) = []; 
    nbrdamagedB(z) = length(XdamagedB);
    percent_damagedB_cells(z)=(length(XdamagedB))/L2(z)*100;
     
    % Set a fixed gate for the population with high FL1-H (SYBRgreen) and 
    % low FL3-H (Propidium iodide) signal. This represent the intact cells.
    % For the dataset used for developing the script, it was found that the
    % line y=0.96*x + 0.6919 was suitable to remove such events from the 
    % log transformed data
    Xintact=X; 
    Xintact(any(Xintact(:,3,:)<=0.96*Xintact(:,4,:) + 0.6919,2),:) = [];
    nbrintact(z) = length(Xintact);
    percent_intact_cells(z)=(length(Xintact))/L2(z)*100; 
    
    % Set a fixed gate for the population of cells that fall inbetween the 
    % damaged B and intact cells, unsing the previous two equations. This 
    % represents the other the damaged cell subpopulation, here called damaged A.
    XdamagedA=X;
    XdamagedA(any(XdamagedA(:,3,:)<=0.908*XdamagedA(:,4,:) + 0.0265,2),:) = [];
    XdamagedA(any(XdamagedA(:,3,:)>=0.96*XdamagedA(:,4,:) + 0.6919,2),:) = [];
    nbrdamagedA(z) = length(XdamagedA);
    percent_damagedA_cells(z)=(length(XdamagedA))/L2(z)*100;

    % Calculate key parameters for the intact cells subpopulation     
    mean_fsc_h_intact(z)= mean(Xintact(:,1));
    mean_ssc_h_intact(z)= mean(Xintact(:,2));
    mode_intact(z)= mode(Xintact(:,1));
    SD_intact (z)= std(Xintact(:,1));
    median_intact(z)= median(Xintact(:,1));
    skewness_intact(z)= mean(Xintact(:,1))-mode(Xintact(:,1));
    rSD1_intact(z)=  mean(Xintact(:,1))-median(Xintact(:,1));
    rSD_intact(z)= abs(mean(Xintact(:,1))-median(Xintact(:,1)))*(1.4826);
    CV_intact(z)= std(Xintact(:,1))/mean(Xintact(:,1));
    rCV_intact(z)= (abs(mean(Xintact(:,1))-median(Xintact(:,1)))*(1.4826))/median(Xintact(:,1));

    %%% Plot the data
    
    % Scatter plot of FSC-H vs. SSC-H for all samples (to visualise light scatter of damagedB and intact cells)
    figure(1);  
    scatter(X(:,1),X(:,2), 2, 'k');
    xlabel('Log FSC-H','FontSize',8); ylabel('Log SSC-H','FontSize',8);
    title(FCMdata(z).fcshrd.filename,'FontSize',8)
    hold on
    scatter(XdamagedB(:,1),XdamagedB(:,2), 2, 'r','filled'); %colored in red
    scatter(Xintact(:,1),Xintact(:,2), 2, 'g','filled'); %colored in green
    scatter(XdamagedA(:,1),XdamagedA(:,2), 2, 'y','filled'); %colored in yellow
    hold off
            
    % Scatter plot of FL1-H vs. FL3-H for all samples (to visualise damagedB and intact cells)         
    figure(2);
    scatter(X(:,3),X(:,4), 3, 'k');
    xlabel('Log FL1-H','FontSize',8); ylabel('Log FL3-H','FontSize',8);
    title(FCMdata(z).fcshrd.filename,'FontSize',8)
    hold on
    scatter(XdamagedB(:,3),XdamagedB(:,4), 1, 'r','filled'); %colored in red
    scatter(Xintact(:,3),Xintact(:,4), 1, 'g','filled'); %colored in green
    scatter(XdamagedA(:,3),XdamagedA(:,4), 2, 'y','filled');%colored in yellow
    hold off
             
end

% concTotal=L3/injection_volume*dilution_factor*1000; %L2 is divided by 20ul, multiplied with 2*1000 (dilution factor) and 1000 (to get per ml instead of per ul) 
% lnconcTotal=log(concTotal);
% concintact=(percent_intact_cells.*L3)/100; %determine number of intact cells in each sample
% lnconcintact=log(concintact);
% concdamagedA=(percent_damagedA_cells.*L2)/100;
% concdamagedB=(percent_damagedB_cells.*L2)/100;
% %timesincestart=hours(timesincestart); %change from duration array to vectors with "hours since start"

concTotal =(L2/injection_volume)*dilution_factor*1000;
concintact = nbrintact/injection_volume*dilution_factor*1000;
concdamagedA = nbrdamagedA/injection_volume*dilution_factor*1000;
concdamagedB = nbrdamagedB/injection_volume*dilution_factor*1000;

%Growth curve, concentration in each gate intact, damagedB, damagedA over time.
figure(5)
time = 0:0.5:(z/2-0.5);
plot(time, concintact,'kx--');
hold on
plot(time, concdamagedA,'bx--');
plot(time, concdamagedB,'rx--');
legend('intact','damagedA','damagedB','Location','northwest');
xlabel('Time (h)');
ylabel('Concentration (cells/ml)');
title('Growth curve');


%Summaries of frequency, cell count in each gate
summary = [percent_intact_cells.', nbrintact.', percent_damagedA_cells.', nbrdamagedA.', percent_damagedB_cells.', nbrdamagedB.',L2.'];
