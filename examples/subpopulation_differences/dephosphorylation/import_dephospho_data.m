% Script for the import of the dephosphorylation data.

filename = 'KM14_KM184-187__ResultsFinal_Cells_final.csv';
delimiter = '\t';
startRow = 2;

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,7,8,9,10,11,12]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
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
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
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
rawNumericColumns = raw(:, [2,7,8,9,10,11,12]);
rawCellColumns = raw(:, [1,3,4,5,6]);


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Allocate imported array to column variable names
PlateID = rawCellColumns(:, 1);
Well = cell2mat(rawNumericColumns(:, 1));
Stain = rawCellColumns(:, 2);
Cond = rawCellColumns(:, 3);
Comp = rawCellColumns(:, 4);
Conc = rawCellColumns(:, 5);
Time = cell2mat(rawNumericColumns(:, 2));
Area = cell2mat(rawNumericColumns(:, 3));
Ch1 = cell2mat(rawNumericColumns(:, 4));
Ch2 = cell2mat(rawNumericColumns(:, 5));
Ch3 = cell2mat(rawNumericColumns(:, 6));
Ch4 = cell2mat(rawNumericColumns(:, 7));


%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;


t = [0,1,4,7,10,13,16,19,22,25,28,31,34,37];
plates = {'KM14_ECM_160808_KM184',...
    'KM14_ECM_160808_KM185',...
    'KM14_ECM_160822_KM186',...
    'KM14_ECM_160822_KM187'};


for r = 1:4
    for i = 1:numel(t)
        ind = find(Time == t(i) & strcmp(Cond,'Ctrl') & strcmp(Comp,'Ctrl') & strcmp(PlateID,plates{r}));
        D.ctrlctrl(r).TrkA{i} = Ch2(ind);
        D.ctrlctrl(r).pErk{i} = Ch3(ind);
    end
end
for r = 1:4
    for i = 1:numel(t)
        ind = find(Time == t(i) & strcmp(Cond,'NGF+GDNF') &strcmp(Comp,'U0126') & strcmp(PlateID,plates{r}));
        D.WFinh(r).TrkA{i} = Ch2(ind);
        D.WFinh(r).pErk{i} = Ch3(ind);
    end
end
for r = 1:4
    for i = 1:numel(t)
        ind = find(Time == t(i) & strcmp(Cond,'Ctrl') & strcmp(Comp,'U0126') & strcmp(PlateID,plates{r}));
        D.ctrlinh(r).TrkA{i} = Ch2(ind);
        D.ctrlinh(r).pErk{i} = Ch3(ind);
    end
end
for r = 1:4
    for i = 1:numel(t)
        ind = find(Time == t(i) & strcmp(Cond,'NGF+GDNF') &strcmp(Comp,'Ctrl') & strcmp(PlateID,plates{r}));
        D.WFctrl(r).TrkA{i} = Ch2(ind);
        D.WFctrl(r).pErk{i} = Ch3(ind);
    end
end
D.t = t;

%save data_dephospho D

