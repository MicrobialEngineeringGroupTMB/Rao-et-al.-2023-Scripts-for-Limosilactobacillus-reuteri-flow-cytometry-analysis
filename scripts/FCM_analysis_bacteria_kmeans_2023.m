%%%%% FCM_analysis_bacteria_kmeans %%%%%

% Written by Julia Thomason 2020-05-31
% Updated by Daniel P. Brink 2023-02-21 for the initial GitHub commit

% Purpose: this script was made to visualise FCM data, and apply k-means 
%          clustering to analyse light scatter data by dynamically setting 
%          gates for each sample based the identified clusters. The script 
%          was developed using data from Limosilactobacillus reuteri 
%          stained with SYBRgreen and Propidium iodide and captured with a 
%          BD Accuri C6+ flow cytometer. While FCS3.1 is a standardized 
%          format, small tweaks might be needed when analysing samples from
%          a different model of flow cytometer. The script is capable of 
%          analysing multiple samples, for example from a time series, by 
%          specifying a folder containing the fcs files.
%
% Calculations: number of cells, frequencies of intact, damagedB and 
%          damagedA cells, concentrations. Table summarizing frequency and
%          events of each gate called "summary". Intact cell population 
%          divided into 5 bin based on pulse width, cell count in each bin 
%          calculated. 
%
% Plots:   light scatter bivariate plot (FSC-H vs SSC-H),
%          viability bivariate plots (FL1-H vs LF3-H),
%          FL1-H Histograms of the intact and damagedA populations,
%          concentration of each population over time.
%          Gated cells are visualised as green (intact), yellow(damagedA), 
%          and red(damagedB). Pulse width histogram, intact population 
%          divided into 5 gates.
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
close all

%Load the fcs files
[FCMdata,number_of_fcs_files_in_folder,channel_labels] = FCSData_to_MData(source_dir);

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

%%Pre-allocate variables for the sake of speed
L1 = zeros(1,number_of_fcs_files_in_folder);
L2 = zeros(1,number_of_fcs_files_in_folder);
L3 = zeros(1,number_of_fcs_files_in_folder);
percent_removed_rows = zeros(1,number_of_fcs_files_in_folder);
sumDistance = zeros(1,number_of_fcs_files_in_folder);
percintact = zeros(1,number_of_fcs_files_in_folder);
percdamagedB = zeros(1,number_of_fcs_files_in_folder);
percdamagedA = zeros(1,number_of_fcs_files_in_folder);
nbrintact = zeros(1,number_of_fcs_files_in_folder);
nbrdamagedB = zeros(1,number_of_fcs_files_in_folder);
nbrdamagedA = zeros(1,number_of_fcs_files_in_folder);
meanFSCintact = zeros(1,number_of_fcs_files_in_folder);
meanSSCintact = zeros(1,number_of_fcs_files_in_folder);
meanwidthintact = zeros(1,number_of_fcs_files_in_folder);

%%% 3. Sample processing %%%
for z=1:number_of_fcs_files_in_folder %for many samples
    
    %% Import the matrix for each sample into a temperary matrix called 'Data'
    clear Data;
    Data = double(FCMdata(z).fcsdat(:,:,1));
    %measure the length of data matrix before cleanup of the raw data
    L1(z) = length(Data); 
    
    %Cleanup of raw data for the sample
    %Remove rows with zeros for FL1-H and FL3-H
    %Remove events below FSC-H 1000, smaller events are not cells
    Data(any(Data(:,9,:)==0,2),:) = [];
    Data(any(Data(:,11,:)==0,2),:) = [];
    Data(any(Data(:,7,:)<=1000,2),:) = [];
    L2(z)=length(Data); %length of data matrix after cleanup (threshold, zeros, etc)
    percent_removed_rows(z)=(L1(z)-L2(z))/L1(z)*100; %how many events are removed from the original file.
    
    %Copy relevant channel number data from the Matrix 'Data' into vectors
    FSC_H=(Data(:,7,:));
    SSC_H=(Data(:,8,:));
    FL1_H=(Data(:,9,:));
    FL3_H=(Data(:,11,:));
    width =(Data(:,13,:));
    X = [log10(FL1_H), log10(FL3_H), log10(FSC_H), log10(SSC_H), width];
    
    %% Use k-means to divide the total population into 3 clusters
    opts = statset('Display','final');
    [idx,C,sumD] = kmeans(X(:,1:2),3,'Distance','cosine','Replicates',20,'Options',opts);
    
    %Total sum of square distance error vector for all clusters
    sumDistance(z) = sum(sumD);
    
    %Assign correct cluster to correct population
    [intactidx, damagedBidx, damagedAidx] = clusterAssign(C);
    Xintact = X(idx==intactidx,:);
    XdamagedB = X(idx==damagedBidx,:);
    XdamagedA = X(idx==damagedAidx,:);
    
    %Helpful for plotting subplots without manual dimension input
    [m,n] = subplots(number_of_fcs_files_in_folder);
    plotVariables = [m,n,z];
    
    %Cleanup of cell debris from intact and damagedA clusters
    [Xintact] = histogramCleanintact(Xintact,plotVariables,FCMdata);
    [XdamagedA] = histogramCleandamagedA(XdamagedA,plotVariables,FCMdata);
    
    %X matrix without cell debris
    X_cells = [XdamagedB; Xintact; XdamagedA];
    
    %Percentages from k-means
    percintact(z) = (length(Xintact)/length(X))*100;
    percdamagedB(z) = (length(XdamagedB)/length(X))*100;
    percdamagedA(z) = (length(XdamagedA)/length(X))*100;
    nbrintact(z) = length(Xintact);
    nbrdamagedB(z) = length(XdamagedB);
    nbrdamagedA(z) = length(XdamagedA);
    L3(z) = nbrintact(z) + nbrdamagedB(z) + nbrdamagedA(z); %number of cells after final clean-up (cell debris)
    meanFSCintact(z) = mean(Xintact(:,3));
    meanSSCintact(z) = mean(Xintact(:,4));
    meanwidthintact(z) = mean(Xintact(:,5));
    
    %% Plot the results
    %Histogram of pulse width divided in 5 gates, returns number of events
    %in each gate
    [LX1,LX2,LX3,LX4,LX5,LXintact] = pulseWidthGate(X_cells,Xintact,plotVariables,FCMdata);
    
    %Scatterplot of FSC vs SSC
    plot_FSC_SSC(X_cells,XdamagedB,Xintact,XdamagedA,plotVariables,FCMdata)
    
    %Scatterplot of FL1 vs FL3 of intact, damagedB, damagedA cells population
    %To include cell debris in figure: replace X_cells with X
    plot_FL1_FL3(X,XdamagedB,Xintact,XdamagedA,plotVariables,FCMdata)

end

%%%                                                        %%%
%%% 4. Cell count estimatation based on the final clusters %%%
%%%                                                        %%%

% The number of events is divided by the injection volume, as specified by 
% the user in section 1,then multiplied with the dilution factor of the 
% sample. Finally, the unit is converted from 1/µl to 1/ml by multiplying 
% with 1000. 

concTotal =(L3/injection_volume)*dilution_factor*1000;
concintact = nbrintact/injection_volume*dilution_factor*1000;
concdamagedA = nbrdamagedA/injection_volume*dilution_factor*1000;
concdamagedB = nbrdamagedB/injection_volume*dilution_factor*1000;

%Summaries of frequency, cell count in each gate
summary = [percintact.', nbrintact.', percdamagedA.', nbrdamagedA.', percdamagedB.', nbrdamagedB.',L2.'];
meanintact = [meanFSCintact; meanSSCintact; meanwidthintact];

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

%table summarizing the number of events in each pulse width gate for all samples
summary_width = [LX1; LX2; LX3; LX4; LX5; LXintact];

%%%                              %%%
%%% 5. Definitions of functions  %%%
%%%                              %%%

function [intactidx, damagedBidx, damagedAidx] = clusterAssign(C)
    %Determination of which cluster is intact, damagedB, damagedA
    %Compare FL1_H centroid locations

    [~, intactidx] = max(C(:,1));
    [~, damagedBidx] = min(C(:,1));
    damagedAidx = find(C(:,1 )== median(C(:,1)));
end

function [m,n] = subplots(number_of_fcs_files_in_folder)
    %Automatic subplot numbering based on the number of files in the source
    %directory
    m = nearest(number_of_fcs_files_in_folder^0.5);
    if number_of_fcs_files_in_folder <= 3
        m = 1;
        n = number_of_fcs_files_in_folder;
    else
        if m*m > number_of_fcs_files_in_folder
            n = m;
        else
            n = m + 1;
        end
    end
end

function [Xintact] = histogramCleanintact(Xintact,plotVariables,FCMdata)
    %Histogram of FL1-H intact or damagedA population to remove additional cell debris
    %If the cleanup-limit is set to higher than 1000 FCS-H this might
    %not be necessary (But other data could be removed accidentally)
    
    m = plotVariables(1);
    n = plotVariables(2);
    z = plotVariables(3);
    figure(1);
    sgtitle('Histogram FL1-H intact');
    subplot(m,n,z);
    h_FL1= histfit(Xintact(:,1,:),100,'kernel'); %histogram with 100 bins fitted with a line
    x_FL1= get(h_FL1(2),'xdata'); %line fitted x vector of histogram bins
    y_FL1= get(h_FL1(2),'ydata'); %line fitted y vector of bin heights
    title(FCMdata(z).fcshrd.filename);
    xlabel('logFL1-H');
    axis([2 6 0 inf]);

    %Local minimum points are found from the histogram.
    %The global minimum is selected & the x-axis value found.
    %If no minimum is found the clean-up limit is set to 0.
    clear y_valleys
    clear y_peak
    min_vector_FL1 = islocalmin(y_FL1);
    max_vector_FL1 = islocalmax(y_FL1,'MaxNumExtrema',1);

    if find(min_vector_FL1)
        %To not remove too much of the intact population only valleys
        %with sufficient separtion are found (half max height)

        y_valleys= y_FL1(min_vector_FL1);
        y_peak = max(y_FL1(max_vector_FL1));
        y_valleys = y_valleys(y_valleys <= y_peak/2);

        %Find minimum valley & remove cell debris below limit
        if find(y_valleys)
            [~,idx_x] = min(y_valleys);
            globalmin_FL1(z) = x_FL1(y_FL1 == y_valleys(idx_x));
        else
            globalmin_FL1(z) = 0;
        end
    else
        globalmin_FL1(z) = 0;
    end

    Xintact(any(Xintact(:,1,:) <= globalmin_FL1(z),2),:) = [];

end

function [XdamagedA] = histogramCleandamagedA(XdamagedA,plotVariables,FCMdata)
    %Histogram of FL1-H intact or damagedA population to remove additional cell debris
    %If the cleanup-limit is set to higher than 1000 FCS-H this might
    %not be necessary (But other data could be removed accidentally)
    m = plotVariables(1);
    n = plotVariables(2);
    z = plotVariables(3);
    figure(2);
    sgtitle('Histogram FL1-H damagedA');
    subplot(m,n,z);
    h_FL1= histfit(XdamagedA(:,1,:),100,'kernel'); %histogram with 100 bins fitted with a line
    x_FL1= get(h_FL1(2),'xdata'); %line fitted x vector of histogram bins
    y_FL1= get(h_FL1(2),'ydata'); %line fitted y vector of bin heights
    title(FCMdata(z).fcshrd.filename);
    xlabel('logFL1-H');
    axis([2 6 0 inf]);
    min_FL1 = islocalmin(y_FL1,'MaxNumExtrema',1);
    if find(min_FL1)
        globalmin_FL1(z) = x_FL1(min_FL1);
    else
        globalmin_FL1(z) = 0;
    end
    XdamagedA(any(XdamagedA(:,1,:) <= globalmin_FL1(z),2),:) = [];
end

function [LX1,LX2,LX3,LX4,LX5,LXintact] = pulseWidthGate(X_cells,Xintact,plotVariables,FCMdata)
    %Divide intact cell poulation into 5 gates based on pulse width
    %histogram minimums. LX1-LX5 are the cell count in each gate.

    m = plotVariables(1);
    n = plotVariables(2);
    z = plotVariables(3);
    figure(6);
    subplot(m,n,z);
    sgtitle('Pulse width - total cell population');
    h_width = histfit(X_cells(:,5),100,'kernel'); %histogram with 100 bins fitted with a line
    x_width = get(h_width(2),'xdata'); %line fitted x vector of histogram bins
    y_width = get(h_width(2),'ydata'); %line fitted y vector of bin heights
    min_vector_width = islocalmin(y_width,'MaxNumExtrema',4); %Finds 4 minimum
    title(FCMdata(z).fcshrd.filename);
    xlabel('Pulse width');
    min_x_width = x_width(min_vector_width); %x-values of min-points
    axis([0 300 0 inf]);
    X1 = Xintact(any(Xintact(:,5,:) <= min_x_width(1),2),:);
    X2 = Xintact(any((Xintact(:,5,:) <= min_x_width(2)) & (Xintact(:,5,:) > min_x_width(1)),2),:);
    X3 = Xintact(any((Xintact(:,5,:) <= min_x_width(3)) & (Xintact(:,5,:) > min_x_width(2)),2),:);
    X4 = Xintact(any((Xintact(:,5,:) <= min_x_width(4)) & (Xintact(:,5,:) > min_x_width(3)),2),:);
    X5 = Xintact(any((Xintact(:,5,:) > min_x_width(4)),2),:);
    LX1(z) = length(X1);
    LX2(z) = length(X2);
    LX3(z) = length(X3);
    LX4(z) = length(X4);
    LX5(z) = length(X5);
    LXintact(z) = length(Xintact);
    figure(7);
    subplot(m,n,z);
    sgtitle('Pulse width - intact cell population');
    histfit(Xintact(:,5),100,'kernel'); %histogram with 100 bins fitted with a line
    xline(min_x_width(1),'-');
    xline(min_x_width(2),'-');
    xline(min_x_width(3),'-');
    xline(min_x_width(4),'-');
    xlabel('Pulse width');
    axis([0 300 0 inf]);
    title(FCMdata(z).fcshrd.filename,'FontSize',8)
    hold off
end

function [] = plot_FSC_SSC(X_cells,XdamagedB,Xintact,XdamagedA,plotVariables,FCMdata)
    %Scatterplot of FSC-H/SSC-H of intact, damagedB, damagedA cells population

    m = plotVariables(1);
    n = plotVariables(2);
    z = plotVariables(3);
    figure(3);
    subplot(m,n,z);
    sgtitle('intact, damagedB, damagedA populations - light scatter bivariate plot');
    scatter(X_cells(:,3),X_cells(:,4), 3, 'ko');
    hold on
    h1 = scatter(XdamagedB(:,3),XdamagedB(:,4),2,'r','filled');
    h2 = scatter(XdamagedA(:,3),XdamagedA(:,4),2,'y','filled');
    h3 = scatter(Xintact(:,3),Xintact(:,4),2,'g','filled');
    xlabel('log(FSC-H)');
    ylabel('log(SSC-H)');
    title(FCMdata(z).fcshrd.filename,'FontSize',8)
    hold off
    end
    
function [] = plot_FL1_FL3(X_cells,XdamagedB,Xintact,XdamagedA,plotVariables,FCMdata)
    %Scatterplot of intact, damagedB, damagedA cells population FL1-H/FL3-H
    %To include cell debris in figure: replace X_cells with X

    m = plotVariables(1);
    n = plotVariables(2);
    z = plotVariables(3);
    figure(4);
    subplot(m,n,z);
    sgtitle('intact, damagedB, damagedA populations - viability bivariate plot after cleanup');
    scatter(X_cells(:,1),X_cells(:,2), 3, 'ko');
    hold on
    h1 = scatter(XdamagedB(:,1),XdamagedB(:,2),2,'r','filled');
    h2 = scatter(XdamagedA(:,1),XdamagedA(:,2),2,'y','filled');
    h3 = scatter(Xintact(:,1),Xintact(:,2),2,'g','filled');
    xlabel('log(FL1-H)');
    ylabel('log(FL3-H)');
    title(FCMdata(z).fcshrd.filename,'FontSize',8)
    hold off
end