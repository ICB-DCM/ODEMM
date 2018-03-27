function ExpC = load_NGF_dose_response_ErkpErk(varargin)
% This file creats the data structure for the 2D dose response data of Erk
% and pErk. For the input/output variables see load_NGF_pErk_kin.m.

%% Import the data
[~, ~, raw] = xlsread('/Users/carolinloos/PhD/ODEMM/ODEMMToolbox/examples/subpopulation_differences/data/KM14_KM101KM103KM105KM107_NGFdoseresponse_ERKpERK_ResultsFinal_Cells.xlsx','Tabelle1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,3,4,5]);
raw = raw(:,[6,7,8,9,10,11,12]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
PlateID = cellVectors(:,1);
Well = cellVectors(:,2);
Stain = cellVectors(:,3);
Cond = cellVectors(:,4);
Comp = cellVectors(:,5);
Conc = data(:,1);
Time = data(:,2);
Area = data(:,3);
Ch1chUCHL1 = data(:,4);
Ch2rbpErk = data(:,5);
Ch3moErk = data(:,6);
Ch4DAPI = data(:,7);

%% Clear temporary variables
clearvars data raw cellVectors R;

repl_names = unique(PlateID); % replicates
u = unique(Conc(~isnan(Conc))); % NGF concentrations, NOTE: for the model 
                                % simulations a scaled concentration of
                                % with 1/25 is used
for k = 1:numel(u)
    j = 1;
    ExpC(k).time = 60;
    ExpC(k).stimulus = u(k);
    for e = 1:numel(repl_names)
        for c = 1:nargin
            str_cond = varargin{c};
            ExpC(k).replicate(j).name{c} = [repl_names{e} ' ' str_cond];
            ExpC(k).replicate(j).measurands = {'Erk','pErk'};
            ind = find(strcmp(PlateID,repl_names{e})&...
                strcmp(Cond,str_cond) &...
                strcmp(Comp, 'NGF')&...
                Conc == u(k));
            ExpC(k).replicate(j).data{c} = [Ch3moErk(ind), Ch2rbpErk(ind)];
            ExpC(k).replicate(j).area{c} = Area(ind);
        end
        j = j+1;
    end
end

