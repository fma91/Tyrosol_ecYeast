function run_ecTyrosol_native()
%RUN_ECTYROSOL_NATIVE Run ecFactory on the Ehrlich-style ecTyrosol_native model.
%
%   Model: ecTyrosol_native.mat (built by ../model/build_ecTyrosol_model_raven.m).
%     ecYeastGEM_batch.mat + native Ehrlich tyrosol branch:
%       4-HPP --ARO10(kcat=1000/s)--> 4-HPAA --ADH7(kcat=1000/s)--> tyrosol
%     No heterologous AAS/AdhE bypass, no chassis deletions.
%
%   Medium: minimal D-glucose ('Min').
%
%   Before ecFactory, model.c is switched to the biomass pseudoreaction (same
%   convention as CellFactory run_predictions.m) so ecFSEOF scans growth flux.
%
%   expYield = 0.49 * WT_yield, WT_yield = 0.48 g/g glucose.
%
%   Outputs (written to ../results/):
%     candidates_L1.txt, candidates_L2.txt, candidates_L3.txt,
%     transporter_targets.txt
%
%   Toolboxes (not modified):
%     ~/Documents/ecFactory/code  (run_ecFactory + GECKO 2.0.3)
%     RAVEN Toolbox + Gurobi on the MATLAB path.

HERE = fileparts(mfilename('fullpath'));
PKG_ROOT = fileparts(HERE);
HOME_DIR = char(java.lang.System.getProperty('user.home'));
ECFACTORY_CODE = fullfile(HOME_DIR, 'Documents', 'ecFactory', 'code');
MODEL_FILE = fullfile(PKG_ROOT, 'model', 'ecTyrosol_native.mat');
RESULTS_FOLDER = fullfile(PKG_ROOT, 'results');
DIARY_FILE = fullfile(PKG_ROOT, 'run_ecTyrosol_native.log');

if exist(DIARY_FILE, 'file'); delete(DIARY_FILE); end
diary(DIARY_FILE); diary on;
diaryCleanup = onCleanup(@() diary('off'));

fprintf('=== ecTyrosol_native ecFactory run ===\n');
fprintf('Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('Model file: %s\n', MODEL_FILE);
fprintf('Results: %s\n', RESULTS_FOLDER);

assert(exist(MODEL_FILE, 'file') == 2, ...
    'Missing model: %s. Run build_ecTyrosol_model_raven.m first.', MODEL_FILE);

geckoLink = fullfile(ECFACTORY_CODE, 'GECKO');
assert(isfolder(geckoLink), ...
    'Missing GECKO at %s. Symlink GECKO 2.0.3 there.', geckoLink);

addpath(genpath(ECFACTORY_CODE));

required = {'solveLP', 'setParam', 'haveFlux'};
missing = required(cellfun(@(f) isempty(which(f)), required));
if ~isempty(missing)
    error('Missing on the MATLAB path: %s.', strjoin(missing, ', '));
end

if ~exist(RESULTS_FOLDER, 'dir'); mkdir(RESULTS_FOLDER); end

WT_YIELD = 0.48;
EXP_YIELD = 0.49 * WT_YIELD;

original_pwd = pwd;
pwdCleanup = onCleanup(@() cd(original_pwd));

raw = load(MODEL_FILE);
fn = fieldnames(raw);
assert(numel(fn) == 1, 'Expected single struct in %s', MODEL_FILE);
ecModel = normalize_model_fields(raw.(fn{1}));

fprintf('  Loaded model: n_rxns=%d, n_mets=%d, n_genes=%d\n', ...
    numel(ecModel.rxns), numel(ecModel.mets), numel(ecModel.genes));

CSname = 'D-glucose exchange (reversible)';
ecModel = changeMedia_batch(ecModel, CSname, 'Min');

if any(strcmpi(ecModel.rxns, 'r_2111'))
    ecModel = setParam(ecModel, 'lb', 'r_2111', 0);
    ecModel = setParam(ecModel, 'ub', 'r_2111', 1000);
end

modelParam = struct();
targetIndex = find(ecModel.c);
assert(numel(targetIndex) == 1, 'Expected exactly one nonzero c entry');
modelParam.rxnTarget = ecModel.rxns{targetIndex};
modelParam.CS_MW = 0.18015;
csIdx = find(strcmpi(ecModel.rxnNames, CSname));
assert(~isempty(csIdx), 'Could not find carbon source rxn "%s"', CSname);
modelParam.CSrxn = ecModel.rxns{csIdx};
grIdx = find(strcmpi(ecModel.rxnNames, 'biomass pseudoreaction'));
assert(~isempty(grIdx), 'Could not find biomass pseudoreaction.');
modelParam.growthRxn = ecModel.rxns{grIdx};

fprintf('  modelParam.rxnTarget = %s (%s)\n', modelParam.rxnTarget, ecModel.rxnNames{targetIndex});
fprintf('  modelParam.CSrxn     = %s\n', modelParam.CSrxn);
fprintf('  modelParam.growthRxn = %s\n', modelParam.growthRxn);

ecModel = setParam(ecModel, 'obj', modelParam.growthRxn, 1);
ecModel = setParam(ecModel, 'lb', modelParam.growthRxn, 0);
ecModel = setParam(ecModel, 'ub', modelParam.growthRxn, 1000);

provModel = setParam(ecModel, 'ub', modelParam.CSrxn, 1);
sol = solveLP(provModel, 1);
if isempty(sol) || ~isfield(sol, 'x') || isempty(sol.x)
    wtBio = NaN;
else
    wtBio = sol.x(find(provModel.c));
end
flux = haveFlux(provModel, 1e-12, modelParam.rxnTarget);
fprintf('  Min sanity: biomass=%g  haveFlux(target)=%d\n', wtBio, flux);
if isnan(wtBio) || wtBio <= 1e-4
    warning('Min medium is infeasible for ecTyrosol_native; ecFactory may fail.');
end

fprintf('  WT_yield=%g, expYield=%g\n', WT_YIELD, EXP_YIELD);

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
%NORMALIZE_MODEL_FIELDS Patch shape/dtype mismatches from scipy.io.savemat.
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
numericFields = {'lb', 'ub', 'c', 'b', 'rev', 'metCharges', 'metComps', ...
    'rxnConfidenceScores', 'MWs'};
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
