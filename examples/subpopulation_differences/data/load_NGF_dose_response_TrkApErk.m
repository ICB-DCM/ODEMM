function ExpC = load_NGF_dose_response_TrkApErk(varargin)
% This file creats the data structure for the 2D dose response data of TrkA
% and pErk. For the input/output variables see load_NGF_pErk_kin.m.

%% Initialize variables.
filename = 'KM14_KM100KM102KM104KM106_TrkApERK_ResultsFinal_Cells.csv';
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
%%

repl_names = {'KM14_ECM_150722_KM100_','KM14_ECM_150722_KM102_','KM14_ECM_150727_KM104_',...
    'KM14_ECM_150727_KM106_'};
% number of doses
n_u = 6;
% number of replicates
n_e = 4;
% considered replicates
I_exp = 1:n_e;
u = [0,0.008,0.04,0.2,1,5];
% label data
% restore data
for k = 1:n_u
    % index set of all measurements for stimulus s
    % initialization
    j = 1;
    % loop over replicates
    ExpC(k).time = 60;
    ExpC(k).stimulus = u(k);
    
    for e = I_exp
        % Additional information
        for c = 1:nargin
         str_cond = varargin{c};
        ExpC(k).replicate(j).name{c} = [repl_names{e} ' ' str_cond];
        ExpC(k).replicate(j).measurands = {'TrkA','pErk'};
        ind = find(strcmp(PlateID,repl_names{e})&...
             strcmp(Cond,str_cond) & strcmp(Comp, 'NGF')&Conc == u(k));
        % fluorescence and size data
        ExpC(k).replicate(j).data{c} = [Ch2(ind), Ch3(ind)];
        ExpC(k).replicate(j).area{c} = Area(ind);
        end
        % update
        j = j+1;
    end
    
end


