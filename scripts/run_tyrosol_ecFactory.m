function run_tyrosol_ecFactory()
%RUN_TYROSOL_ECFACTORY Genome-scale strain design for tyrosol production.
%
%   Pipeline:
%     1. Load ecTyrosol.mat (built by ../model/build_ecTyrosol_model_raven.m).
%     2. Apply minimal medium with D-glucose as carbon source.
%     3. Run ecFactory (GECKO 2.x + RAVEN) to predict gene targets at three
%        filtering levels (L1, L2, L3) plus transporter reactions.
%
%   Model assumptions: docs/METHODS.md and model/build_ecTyrosol_model_raven.m.
%   This script does not modify any toolbox code.
%
%   Required on the MATLAB path:
%     - ~/Documents/ecFactory/code  (GECKO 2.0.3 via code/GECKO symlink)
%     - RAVEN Toolbox, Gurobi
%
%   Outputs (written to ../results/):
%     candidates_L1.txt, candidates_L2.txt, candidates_L3.txt,
%     transporter_targets.txt

HERE = fileparts(mfilename('fullpath'));
PKG_ROOT = fileparts(HERE);
HOME_DIR = char(java.lang.System.getProperty('user.home'));
ECFACTORY_CODE = fullfile(HOME_DIR, 'Documents', 'ecFactory', 'code');
MODEL_FILE = fullfile(PKG_ROOT, 'model', 'ecTyrosol.mat');
RESULTS_FOLDER = fullfile(PKG_ROOT, 'results');
DIARY_FILE = fullfile(PKG_ROOT, 'run_tyrosol_ecFactory.log');

if exist(DIARY_FILE, 'file'); delete(DIARY_FILE); end
diary(DIARY_FILE); diary on;
diaryCleanup = onCleanup(@() diary('off'));

fprintf('=== Tyrosol strain design (ecFactory) ===\n');
fprintf('Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('Model: %s\n', MODEL_FILE);
fprintf('Results: %s\n', RESULTS_FOLDER);

assert(exist(MODEL_FILE, 'file') == 2, ...
    'Missing %s. Run build_ecTyrosol_model_raven.m first.', MODEL_FILE);
assert(isfolder(fullfile(ECFACTORY_CODE, 'GECKO')), ...
    'GECKO 2.0.3 required at %s/GECKO', ECFACTORY_CODE);

addpath(genpath(ECFACTORY_CODE));
required = {'solveLP', 'setParam', 'haveFlux'};
missing = required(cellfun(@(f) isempty(which(f)), required));
if ~isempty(missing)
    error('Missing on path: %s', strjoin(missing, ', '));
end
if ~exist(RESULTS_FOLDER, 'dir'); mkdir(RESULTS_FOLDER); end

% ecFactory yield scan parameters (CellFactory-ecYeastGEM convention)
WT_YIELD = 0.48;          % g biomass / g glucose on minimal medium
EXP_YIELD = 0.49 * WT_YIELD;

original_pwd = pwd;
pwdCleanup = onCleanup(@() cd(original_pwd));

raw = load(MODEL_FILE);
fn = fieldnames(raw);
assert(numel(fn) == 1, 'Expected one struct in %s', MODEL_FILE);
ecModel = normalize_model_fields(raw.(fn{1}));

fprintf('  Model size: %d reactions, %d metabolites, %d genes\n', ...
    numel(ecModel.rxns), numel(ecModel.mets), numel(ecModel.genes));

% Medium: minimal glucose (assumption — de novo aromatic biosynthesis active)
CSname = 'D-glucose exchange (reversible)';
ecModel = changeMedia_batch(ecModel, CSname, 'Min');

% Relax optional minimum-growth constraint on r_2111 if present
if any(strcmpi(ecModel.rxns, 'r_2111'))
    ecModel = setParam(ecModel, 'lb', 'r_2111', 0);
    ecModel = setParam(ecModel, 'ub', 'r_2111', 1000);
end

modelParam = struct();
targetIndex = find(ecModel.c);
assert(numel(targetIndex) == 1, 'Product objective must be unique in model.c');
modelParam.rxnTarget = ecModel.rxns{targetIndex};
modelParam.CS_MW = 0.18015;
modelParam.CSrxn = ecModel.rxns{strcmpi(ecModel.rxnNames, CSname)};
modelParam.growthRxn = ecModel.rxns{strcmpi(ecModel.rxnNames, 'biomass pseudoreaction')};

fprintf('  Product reaction : %s (%s)\n', modelParam.rxnTarget, ecModel.rxnNames{targetIndex});
fprintf('  Carbon source    : %s\n', modelParam.CSrxn);
fprintf('  Growth reaction  : %s\n', modelParam.growthRxn);

% ecFSEOF requires biomass as the active objective during the scan
ecModel = setParam(ecModel, 'obj', modelParam.growthRxn, 1);
ecModel = setParam(ecModel, 'lb', modelParam.growthRxn, 0);
ecModel = setParam(ecModel, 'ub', modelParam.growthRxn, 1000);

provModel = setParam(ecModel, 'ub', modelParam.CSrxn, 1);
sol = solveLP(provModel, 1);
if isempty(sol) || ~isfield(sol, 'x') || isempty(sol.x)
    fprintf('  Warning: LP check returned empty solution on minimal medium.\n');
else
    wtBio = sol.x(find(provModel.c));
    flux = haveFlux(provModel, 1e-12, modelParam.rxnTarget);
    fprintf('  LP check: max biomass = %g, product flux feasible = %d\n', wtBio, flux);
end

fprintf('  expYield = %.4f (0.49 x WT_yield %.2f)\n', EXP_YIELD, WT_YIELD);

ecfseof_results_dir = fullfile(ECFACTORY_CODE, 'GECKO', 'geckomat', 'utilities', 'ecFSEOF', 'results');
if exist(ecfseof_results_dir, 'dir')
    try; rmdir(ecfseof_results_dir, 's'); catch; end
end

cd(ECFACTORY_CODE);
try
    [~, candidates, step] = run_ecFactory(ecModel, modelParam, EXP_YIELD, RESULTS_FOLDER, false);
    cd(original_pwd);
    fprintf('\n  ecFactory finished at step %d (%d candidates).\n', step, height(candidates));
catch ME
    cd(original_pwd);
    rethrow(ME);
end

end


function model = normalize_model_fields(model)
%NORMALIZE_MODEL_FIELDS Convert Python-exported .mat fields for RAVEN/Gurobi.
matrixFields = {'S', 'rxnGeneMat'};
fn = fieldnames(model);
for k = 1:numel(fn)
    f = fn{k};
    if ismember(f, matrixFields); continue; end
    v = model.(f);
    sz = size(v);
    if numel(sz) ~= 2; continue; end
    if sz(1) == 1 && sz(2) > 1
        if ischar(v) || isstring(v)
            v = cellstr(v);
        end
        model.(f) = v(:);
    elseif ischar(v) && sz(1) > 1
        model.(f) = cellstr(v);
    end
end
numericFields = {'lb', 'ub', 'c', 'b', 'rev', 'metCharges', 'metComps', 'rxnConfidenceScores', 'MWs'};
for i = 1:numel(numericFields)
    f = numericFields{i};
    if isfield(model, f)
        model.(f) = full(double(model.(f)(:)));
    end
end
if isfield(model, 'S') && ~issparse(model.S)
    model.S = sparse(model.S);
end
if isfield(model, 'rxnGeneMat') && ~issparse(model.rxnGeneMat)
    model.rxnGeneMat = sparse(model.rxnGeneMat);
end
end
