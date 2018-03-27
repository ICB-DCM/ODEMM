function ExpC = load_NGF_pErk_kin(varargin)
% This file creats the data structure for the 1D NGF-kinetics of pErk.
%
% USAGE:
% ExpC = load_NGF_pErk_kin('PDL') \n
% ExpC = load_NGF_pErk_kin('PDL', 'ColI')
%
% Parameters:
% varargin: strings indicating the extracellular scaffold, e.g. 'PDL' or
% 'ColI'
%
% Return values:
% ExpC:  struct with fields
%   * name: string specifying the replicateal conditions
%   * label: label used in plots for ticks
%   * labelname: label used in plots as label
%   * replicate(j): struct with fields
%      - name: {'string specifying the 1st replicate/condition','string
%      specifying the 2nd replicate/condition',...}
%       - measurands: {name of measurand 1,name of measurand 2,
%                            ...,name of measurand m}
%       - data:  {n_D1 x m matrix under condition 1, n_D2 x m matrix under
%          condition 2,...,n_Dnc x m matrix under condition n_c}
%         (One row represents one observed cell with the
%          data in the order of the measurands. The different
%          rows provide measurement data for different cells.)

%% Import the data
[~, ~, raw] = xlsread('KM14_KM28KM31KM34_NGFkinetic_pErk_ResultsFinal_Cells.xlsx','Tabelle1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,3,4,5,6]);
raw = raw(:,[7,8,9,10,11]);

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
Conc = cellVectors(:,6);
Time = data(:,1);
Area = data(:,2);
Ch2rbpErk = data(:,4);

%% Clear temporary variables
clearvars data raw cellVectors R;

repl_names = unique(PlateID); % replicate names
t = [0,1,5,15,30,60,120]; % time points

for k = 1:numel(t)
    j = 1;
    ExpC(k).time = t(k);
    ExpC(k).stimulus = 1;
    for e = 1:numel(repl_names)
        ExpC(k).replicate(j).name = repl_names{e};
        ExpC(k).replicate(j).measurands = {'pErk'};
        
        for c = 1:nargin
            str_cond = varargin{c};
            ExpC(k).replicate(j).name_condition{c} = str_cond;
            ExpC(k).replicate(j).data{c} = Ch2rbpErk(strcmp(PlateID,repl_names{e}) &...
                 strcmp(Cond,str_cond) &...
                 strcmp(Comp, 'NGF') &...
                 Time == t(k));
        end
        j = j+1;
    end
end


