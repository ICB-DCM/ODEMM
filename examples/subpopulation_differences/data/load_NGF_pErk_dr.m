function ExpC = load_NGF_pErk_dr(varargin)
% This file creats the data structure for the 1D NGF-dose response of pErk.
% For the input/output parameters see load_NGF_pErk_kin.m.

%% Import the data
[~, ~, raw] = xlsread('KM14_KM60KM62KM67KM68_NGFdoseresponse_pErk_ResultsFinal_Cells.xlsx','Tabelle1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,3,4,5]);
raw = raw(:,[6,7,8,9,10,11]);

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
Ch3DAPI = data(:,6);

%% Clear temporary variables
clearvars data raw cellVectors R;

repl_names = unique(PlateID); % replicate names
u = unique(Conc(~isnan(Conc))); % NGF concentrations, NOTE: for the model 
                                % simulations a scaled concentration of
                                % with 1/25 is used

for k = 1:numel(u)
    j = 1;
    ExpC(k).time = 60;
    ExpC(k).stimulus = u(k);
    for e = 1:numel(repl_names)
        ExpC(k).replicate(j).name = repl_names{e};
        ExpC(k).replicate(j).measurands = {'pErk'};
        for c= 1:nargin
            str_cond = varargin{c};
            ExpC(k).replicate(j).name_condition{c} = str_cond;
            ExpC(k).replicate(j).data{c} = Ch2rbpErk(strcmp(PlateID,repl_names{e})&...
                 strcmp(Cond,str_cond) &...
                 strcmp(Comp, 'NGF')&...
                 Conc == u(k));
        end
        j = j+1;
    end
end


