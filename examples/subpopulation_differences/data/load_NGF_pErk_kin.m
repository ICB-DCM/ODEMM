function ExpC = load_NGF_pErk_kin(varargin)
% This file creats the data structure for the 1D kinetic data.
%
% USAGE:
% ExpC = load_NGF_pErk_kin(str1) \n
% ExpC = load_NGF_pErk_kin(str1,str2)
%
% Parameters:
% varargin:
% * strings indicating the extracellular scaffold used ("LYS","Col
%
% Return values:
% ExpC(i)
%   .name = 'string specifying the replicateal consitions'
%   .label = 'label used in plots for ticks'
%   .labelname = 'label used in plots as label'
%   .replicate(j)
%      .name = {'string specifying the 1st replicate/condition','string
%      specifying the 2nd replicate/condition',...}
%      .measurands = {'name of measurand 1','name of measurand 2',
%                            ...,'name of measurand m'}
%      .data =  {n_D1 x m matrix under condition 1, n_D2 x m matrix under
%          condition 2,...,n_Dnc x m matrix under condition n_c}
%         (One row represents one observed cell with the
%          data in the order of the measurands. The different
%          rows provide measurement data for different cells.)
%

%% Initialize variables.
filename = 'KM14_KM28_31_34_ResultsFinal_Cells.csv';
delimiter = '\t';
startRow = 2;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3,6,7,8,9,10,11,12]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [2,3,6,7,8,9,10,11,12]);
rawCellColumns = raw(:, [1,4,5]);


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Allocate imported array to column variable names
PlateID = rawCellColumns(:, 1);
Well = cell2mat(rawNumericColumns(:, 1));
Stain = cell2mat(rawNumericColumns(:, 2));
Cond = rawCellColumns(:, 2);
Comp = rawCellColumns(:, 3);
Conc = cell2mat(rawNumericColumns(:, 3));
Time = cell2mat(rawNumericColumns(:, 4));
Area = cell2mat(rawNumericColumns(:, 5));
Ch1 = cell2mat(rawNumericColumns(:, 6));
Ch2 = cell2mat(rawNumericColumns(:, 7));
Ch3 = cell2mat(rawNumericColumns(:, 8));
Ch4 = cell2mat(rawNumericColumns(:, 9));


%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;

repl_names = {'KM14_ECM_140901_KM28_','KM14_ECM_140903_KM31_','KM14_ECM_140915_KM34_'};
% number of time points
n_t = 7;
% number of replicates
n_e = 3;
% considered experiments
I_exp = 1:n_e;
t = [0,1,5,15,30,60,120];
% label data
% restore data
for k = 1:n_t
    % index set of all measurements for stimulus s
    % initialization
    j = 1;
    % loop over experiments
    ExpC(k).time = t(k);
    ExpC(k).stimulus = 1;
    for e = I_exp
        % Additional information
        ExpC(k).replicate(j).name = repl_names{e};
        ExpC(k).replicate(j).measurands = {'pErk'};
        
        for c = 1
            str_cond = 'LYS';
            ExpC(k).replicate(j).name_condition{c} = str_cond;
            ind = find(strcmp(PlateID,repl_names{e})&...
                 strcmp(Cond,str_cond) & strcmp(Comp, 'NGF')&Time == t(k));
             numel(ind)
            ExpC(k).replicate(j).data{c} = Ch2(ind);
            ExpC(k).replicate(j).area{c} = Area(ind);
        end
        % update
        j = j+1;
    end
end
