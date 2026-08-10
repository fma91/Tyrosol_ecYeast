function run_ecTyrosol_native_no_mito()
%RUN_ECTYROSOL_NATIVE_NO_MITO ecFactory (FSEOF + EUVA) on the rho0 / no-mito model.
%
%   Same algorithm and settings as run_ecTyrosol_native.m, with one difference:
%   the seven mtDNA-encoded genes are knocked out before the scan (same
%   noMito helper used by tyrosol_envelopes.m and tyrosol_dfba.m).
%
%   Model: model/ecTyrosol_native.mat (or ecTyrosol.mat alias)
%   Medium: Min, D-glucose
%   Yield scan: WT_yield = max biomass yield of the no-mito model on Min;
%     expYield = 0.49 * WT_yield (same CellFactory convention as the WT run).
%     Do NOT reuse WT_yield=0.48 — that forces an infeasible biomass floor.
%
%   Outputs (does NOT overwrite wild-type results/):
%     results/no_mito/candidates_L1.txt
%     results/no_mito/candidates_L2.txt
%     results/no_mito/candidates_L3.txt
%     results/no_mito/transporter_targets.txt
%
%   Usage:
%     addpath('scripts')
%     run_ecTyrosol_native_no_mito

HERE = fileparts(mfilename('fullpath'));
PKG_ROOT = fileparts(HERE);
HOME_DIR = char(java.lang.System.getProperty('user.home'));
ECFACTORY_CODE = fullfile(HOME_DIR, 'Documents', 'ecFactory', 'code');
MODEL_FILE = fullfile(PKG_ROOT, 'model', 'ecTyrosol_native.mat');
if exist(MODEL_FILE, 'file') ~= 2
    MODEL_FILE = fullfile(PKG_ROOT, 'model', 'ecTyrosol.mat');
end
RESULTS_FOLDER = fullfile(PKG_ROOT, 'results', 'no_mito');
DIARY_FILE = fullfile(PKG_ROOT, 'run_ecTyrosol_native_no_mito.log');

if exist(DIARY_FILE, 'file'); delete(DIARY_FILE); end
diary(DIARY_FILE); diary on;
diaryCleanup = onCleanup(@() diary('off'));

fprintf('=== ecTyrosol_native ecFactory run (no-mito / rho0) ===\n');
fprintf('Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('Model file: %s\n', MODEL_FILE);
fprintf('Results: %s\n', RESULTS_FOLDER);

assert(exist(MODEL_FILE, 'file') == 2, ...
    'Missing model: %s. Run build_ecTyrosol_model_raven.m first.', MODEL_FILE);

geckoLink = fullfile(ECFACTORY_CODE, 'GECKO');
assert(isfolder(geckoLink), ...
    'Missing GECKO at %s. Symlink GECKO 2.0.3 there.', geckoLink);

addpath(genpath(ECFACTORY_CODE));

required = {'solveLP', 'setParam', 'haveFlux', 'removeGenes'};
missing = required(cellfun(@(f) isempty(which(f)), required));
if ~isempty(missing)
    error('Missing on the MATLAB path: %s.', strjoin(missing, ', '));
end

if ~exist(RESULTS_FOLDER, 'dir'); mkdir(RESULTS_FOLDER); end

% Yield scan: same CellFactory convention (expYield = 0.49 * WT_yield), but
% WT_yield must be the no-mito max biomass yield on Min. The wild-type value
% 0.48 forces V_bio ~0.042 h-1 in EUVA, which is infeasible on rho0 (max mu
% ~0.012 h-1 at unit glucose uptake).
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

% rho0 / no-mito background (identical gene list to envelopes and dFBA)
ecModel = noMito(ecModel);
fprintf('  Applied no-mito: KO of Q0045,Q0080,Q0085,Q0105,Q0130,Q0250,Q0275\n');
fprintf('  After no-mito: n_rxns=%d, n_mets=%d, n_genes=%d\n', ...
    numel(ecModel.rxns), numel(ecModel.mets), numel(ecModel.genes));

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
fprintf('  Min + no-mito sanity: biomass=%g  haveFlux(target)=%d\n', wtBio, flux);
if isnan(wtBio) || wtBio <= 1e-4
    error('Min + no-mito is infeasible; cannot calibrate WT_yield for ecFactory.');
end

% g biomass / g glucose at unit glucose uptake (mmol/gDW/h)
WT_YIELD = wtBio / modelParam.CS_MW;
EXP_YIELD = 0.49 * WT_YIELD;
fprintf('  Calibrated WT_yield=%.4f g/g (no-mito max), expYield=%.4f\n', ...
    WT_YIELD, EXP_YIELD);
fprintf('  (Wild-type run used WT_yield=0.48; that value is infeasible here.)\n');

ecfseof_results_dir = fullfile(ECFACTORY_CODE, 'GECKO', 'geckomat', 'utilities', 'ecFSEOF', 'results');
if exist(ecfseof_results_dir, 'dir')
    try; rmdir(ecfseof_results_dir, 's'); catch; end
end

cd(ECFACTORY_CODE);
try
    [~, candidates, step] = run_ecFactory(ecModel, modelParam, EXP_YIELD, RESULTS_FOLDER, false);
    cd(original_pwd);
    fprintf('\n  ecFactory (no-mito) finished at step %d (%d candidates).\n', step, height(candidates));
    fprintf('  Tables written to %s\n', RESULTS_FOLDER);
catch ME
    cd(original_pwd);
    rethrow(ME);
end

end


function m = noMito(m)
%NOMITO Knock out the seven mtDNA-encoded genes (rho0 / petite background).
for g = {'Q0045', 'Q0080', 'Q0085', 'Q0105', 'Q0130', 'Q0250', 'Q0275'}
    if any(strcmpi(m.genes, g{1})) || ...
            (isfield(m, 'enzGenes') && any(strcmpi(m.enzGenes, g{1})))
        m = removeGenes(m, g{1}, false, false, false);
    end
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
