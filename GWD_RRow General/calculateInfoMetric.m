function InfoVal = calculateInfoMetric(inputData,method,runUnique)

% Takes in a single vector (N x 1) or pair of vectors (N x 2) for
% calculating Entropy or the Mutual Information. Input a str of 'Entropy'
% or 'PairedMI' or 'JointMI#' for method to identify which calculation to make.


% GWD Sept 2020 Adapted/Simplified from the Neuroscience Information
% Toolbox (instinfo.m) written by Nick Timme:
%
% Email: nicholas.m.timme@gmail.com
% August 2014; Last revision: 12-May-2016
% https://github.com/nmtimme/Neuroscience-Information-Theory-Toolbox

% Currently ported functionalitly: Single variable entropy, Paired Mutual
% Information, Joint Mutual Information (2 groups, unlimited entries in
% each group)
%
% NOTE: Significance testing via MonteCarlo simulations is not currently
% supported internally, shuffle your data outside and then pump it through
% here a bunch.

% Default to running unique to ensure that your data are sequential
% integers starting from 1. However, as this takes the most time of the
% analysis you can run the beforehand on the data and bypass this option.


if contains(method,'JointMI')
    transIdx = str2double(strrep(method,'JointMI',''));
    method = 'JointMI';
    assert(size(inputData,2) > transIdx,'You dont have enough data elements to split as requested')
end

switch method
    case 'Entropy'
        
        % Remove any nan entries from the input data
        x = inputData(~isnan(inputData));
        
        % There must be at least 2 entries to compute Information metrics
        if length(x) < 2
            InfoVal = nan;
            return
        end
        
        % Ensure data run up from 1 as integers
        if ~exist('runUnique','var') || runUnique 
            [~,~,x] = unique(x);
        end
        
        % Make a counts matrix
        Counts = accumarray({x},ones(size(x)),[max(x),1]);
        
        % Calculate the entropy
        InfoVal = EntropyY(Counts);
        
    case 'JointEntropy'
        
        error('Joint Entropy not yet implemented')
        % % %
        % % %         [B,I,x] = unique(x,'rows');
        % % %
        % % %     % Error check the integer states
        % % %     if ~isequal(unique(x),(1:length(unique(x))))
        % % %         [B,I,x] = unique(x);
        % % %     end
        % % %
        % % %     % Make a counts matrix
        % % %     Counts = accumarray({x},ones(size(x)),[length(unique(x)),1]);
        % % %
        % % %     % Calculate the entropy
        % % %     InfoVal = EntropyY(Counts);
        % % %
    case 'PairedMI'
        
        x = inputData(:,1);
        y = inputData(:,2);
        
        % Remove nans
        noData = any(isnan(inputData),2);
        x(noData) = [];
        y(noData) = [];
        
        % There must be at least 2 entries to compute Information metrics
        if length(x) < 2
            InfoVal = nan;
            return
        end
        
        % Ensure data run up from 1 as integers
        if ~exist('runUnique','var') || runUnique
            [~,~,x] = unique(x);
            [~,~,y] = unique(y);
        end
        
        % Make a counts matrix
        Counts = accumarray({x,y},ones(size(x)),[max(x),max(y)]);
        
        % Calculate the pairwise mutual information
        InfoVal = MutualInfo(Counts);
        
    case 'JointMI'
        
        x = inputData(:,1:transIdx);
        y = inputData(:,transIdx+1:end);
        
        % Remove any entries that have nans
        noData = any(isnan(inputData),2);
        x(noData,:) = [];
        y(noData,:) = [];
        
        % There must be at least 2 entries to compute Information metrics
        if length(x) < 2
            InfoVal = nan;
            return
        end
        
        % Convert the data into integers
        if ~exist('runUnique','var') || runUnique
            [~,~,x] = unique(x,'rows');
            [~,~,y] = unique(y,'rows');
        end
        
        % Make a counts matrix
        Counts = accumarray({x,y},ones(size(x)),[max(x),max(y)]);
        
        % Calculate the pairwise mutual information
        InfoVal = MutualInfo(Counts);
        
    case 'TransferEntropy'
        %%% CURRENTLY TESTING %%%
        %%% CURRENTLY TESTING %%%
        
        y_Future = inputData(:,1); % Future/Predicted Y Data
        y_Past = inputData(:,2); % Previous/Reference Y Data
        x = inputData(:,3); % Alternative/External/Other predictor Data
        
        % Remove nans
        noData = any(isnan(inputData),2);
        x(noData) = [];
        y_Past(noData) = [];
        y_Future(noData) = [];
        
        % There must be at least 2 entries to compute Information metrics
        if length(x) < 2
            InfoVal = nan;
            return
        end
        
        % Ensure data run up from 1 as integers
        if ~exist('runUnique','var') || runUnique
            [~,~,x] = unique(x);
            [~,~,y_Past] = unique(y_Past);
            [~,~,y_Future] = unique(y_Future);
        end
        
        % Make a counts matrix
        Counts = accumarray({y_Future,y_Past,x},ones(size(x)),[max(y_Future),max(y_Past),max(x)]);
        
        % Calculate the pairwise mutual information
        InfoVal = TE2(Counts);
        
        %%% CURRENTLY TESTING %%%
        %%% CURRENTLY TESTING %%%
end