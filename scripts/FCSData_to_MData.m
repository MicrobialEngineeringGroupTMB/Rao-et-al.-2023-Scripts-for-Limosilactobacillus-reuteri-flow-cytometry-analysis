function [FCMdata,number_of_fcs_files_in_folder, channel_labels] = FCSData_to_MData(source_dir)
%%%% Converting FCS files into MAT %%%%

%%% Wrapper script to load fcs files using the fca_readfcs function
%
% Dependencies: fca_readfcs by Laszlo Balkay, available from from the Matlab Central
% Full reference:
% Laszlo Balkay (2023). fca_readfcs (https://www.mathworks.com/matlabcentral/fileexchange/9608-fca_readfcs), 
% MATLAB Central File Exchange. Retrieved January 27, 2023. 

%Save the start time of the run to variable
tstart=tic;

%Load all fcs files in source_dir
source_files = dir(fullfile(source_dir,'*.fcs'));
number_of_fcs_files_in_folder=length(source_files);

%loop over all fcs files in dir, load them using fca_readfcs 
% and convert the data to a Matlab data structure
for z=1:number_of_fcs_files_in_folder
    [fcsdat, fcshdr, fcsdatscaled]=fca_readfcs(fullfile(source_dir,source_files(z).name));
    FCMdata(z).fcsdat(:,:,1)=fcsdat;
    FCMdata(z).fcshrd(:,:,1)=fcshdr;
    FCMdata(z).datscaled(:,:,1)=fcsdatscaled;
           
end

% Get the channel labele
channel_labels={fcshdr.par.name}

%Calculate the elapsed time
tElapsedtime=toc(tstart);

%Print a summary of the number of analyzed fcs files and elapsed time 
table(number_of_fcs_files_in_folder, tElapsedtime)
