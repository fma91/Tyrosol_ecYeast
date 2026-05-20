function build_ecTyrosol_model_raven(outFile, baseFile)
%BUILD_ECTYROSOL_MODEL_RAVEN Build ecTyrosol.mat with RAVEN Toolbox.
%
%   Reproduces the enzyme-constrained tyrosol model from ecYeastGEM_batch.mat
%   using standard RAVEN model-editing functions (addMets, addRxns, setParam).
%   Run scripts/run_tyrosol_ecFactory.m next to reproduce strain-design targets
%   with ecFactory (GECKO 2.0.3 + RAVEN + Gurobi).
%
%   To adapt this workflow to another product, edit metsToAdd and rxnsToAdd
%   below (metabolite IDs, stoichiometry, bounds, grRules, objective reaction).
%   Enzyme-constrained steps use the GECKO pattern: protein pseudo-metabolite
%   in S with coefficient 1 / (kcat * 3600); here kcat = 1000 1/s (assumption A4).
%
%   Modeling assumptions A1–A6: docs/METHODS.md
%
%   Usage:
%     cd model
%     build_ecTyrosol_model_raven
%
%   Environment variable (optional): ECYEASTGEM_BATCH — path to base .mat file
%
%   Requires RAVEN Toolbox (addMets, addRxns, setParam) on the MATLAB path.

if nargin < 1 || isempty(outFile)
    outFile = fullfile(fileparts(mfilename('fullpath')), 'ecTyrosol.mat');
end
if nargin < 2 || isempty(baseFile)
    baseFile = getenv('ECYEASTGEM_BATCH');
    if isempty(baseFile)
        homeDir = char(java.lang.System.getProperty('user.home'));
        baseFile = fullfile(homeDir, 'Documents', 'CellFactory-ecYeastGEM', ...
            'ModelFiles', 'ecYeastGEM_batch.mat');
    end
end

assert(exist(baseFile, 'file') == 2, 'Base model not found: %s', baseFile);
assert(~isempty(which('addRxns')), ...
    'RAVEN addRxns not on path. Add RAVEN Toolbox before running this script.');
assert(~isempty(which('addMets')), ...
    'RAVEN addMets not on path. Add RAVEN Toolbox before running this script.');

% Assumption A4
kcat = 1000;
sAro10 = 1 / (kcat * 3600);
sAdh7 = 1 / (kcat * 3600);

fprintf('=== build_ecTyrosol_model_raven ===\n');
fprintf('Base : %s\n', baseFile);
fprintf('Out  : %s\n', outFile);

raw = load(baseFile);
fn = fieldnames(raw);
assert(numel(fn) == 1, 'Expected one model struct in %s', baseFile);
varName = fn{1};
model = raw.(varName);

fprintf('Loaded %s (%d rxns, %d mets)\n', varName, numel(model.rxns), numel(model.mets));

% Clear any previous product objective before extending (matches Python build)
model.c(:) = 0;

% --- Assumption A6: new metabolites (addMets) ---
metsToAdd.mets = {'s_4hpaa_c', 's_tyrosol_c', 's_tyrosol_e'};
metsToAdd.metNames = {'4-hydroxyphenylacetaldehyde', 'tyrosol', 'tyrosol'};
metsToAdd.compartments = {'c', 'c', 'e'};
if isfield(model, 'metFormulas')
    metsToAdd.metFormulas = {'C8H8O2', 'C8H10O2', 'C8H10O2'};
end
model = addMets(model, metsToAdd);

% Existing metabolites used in new reactions
mH = metId(model, 's_0794');
mCO2 = metId(model, 's_0456');
mNADPH = metId(model, 's_1212');
mNADP = metId(model, 's_1207');
mHpp = metId(model, findMetByNameComp(model, '3-(4-hydroxyphenyl)pyruvate', 'cytoplasm'));
mProtAro10 = metId(model, 'prot_Q06408');  % ARO10, assumption A2
mProtAdh7 = metId(model, 'prot_P25377');   % ADH7, assumption A3

% --- Assumptions A2–A6: four reactions (addRxns) ---
rxnsToAdd.rxns = {'new_aro10_HPP', 'new_adh7_tyrosol', 'new_tyrosol_t', 'new_tyrosol_ex'};
rxnsToAdd.rxnNames = { ...
    '4-hydroxyphenylpyruvate decarboxylase (ARO10)', ...
    '4-hydroxyphenylacetaldehyde reductase (ADH7)', ...
    'tyrosol transport', ...
    'tyrosol exchange'};
rxnsToAdd.mets = { ...
    {mH, mHpp, mProtAro10, 's_4hpaa_c', mCO2}, ...
    {mH, mNADPH, 's_4hpaa_c', mProtAdh7, mNADP, 's_tyrosol_c'}, ...
    {'s_tyrosol_c', 's_tyrosol_e'}, ...
    {'s_tyrosol_e'}};
rxnsToAdd.stoichCoeffs = { ...
    [-1, -1, -sAro10, 1, 1], ...
    [-1, -1, -1, -sAdh7, 1, 1], ...
    [-1, 1], ...
    [-1]};
rxnsToAdd.lb = [0; 0; 0; 0];
rxnsToAdd.ub = [1000; 1000; 1000; 1000];
rxnsToAdd.grRules = {'YDR380W'; 'YCR105W'; ''; ''};
if isfield(model, 'eccodes')
    rxnsToAdd.eccodes = {'4.1.1.43'; '1.1.1.90'; ''; ''};
end
if isfield(model, 'subSystems')
    rxnsToAdd.subSystems = { ...
        'sce00350  Tyrosine metabolism', ...
        'sce00350  Tyrosine metabolism', ...
        'sce04147  Exosome', ...
        ''};
end
if isfield(model, 'rxnConfidenceScores')
    rxnsToAdd.rxnConfidenceScores = [2; 2; 2; 2];
end

model = addRxns(model, rxnsToAdd, 1);

% Assumption A6: product exchange is the sole objective
if exist('setParam', 'file')
    model = setParam(model, 'obj', 'new_tyrosol_ex', 1);
else
    model.c(:) = 0;
    model.c(strcmp(model.rxns, 'new_tyrosol_ex')) = 1;
end

if isfield(model, 'description')
    model.description = [ ...
        'ecTyrosol: ecYeastGEM_batch extended for tyrosol production via the ', ...
        'Ehrlich pathway (ARO10, ADH7) with enzyme-constrained arms (kcat=1000/s). ', ...
        'Built with RAVEN addMets/addRxns.'];
end

[outDir, ~, ~] = fileparts(outFile);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end
save(outFile, varName, 'model', '-v7');
fprintf('Saved %s (%d rxns, %d mets)\n', outFile, numel(model.rxns), numel(model.mets));
fprintf('Product objective: new_tyrosol_ex\n');
fprintf('Next: addpath(''../scripts''); run_tyrosol_ecFactory\n');

end

% -------------------------------------------------------------------------
function id = metId(model, idxOrId)
if ischar(idxOrId) || isstring(idxOrId)
    id = char(idxOrId);
    return;
end
id = model.mets{idxOrId};
end

function idx = findMetByNameComp(model, metName, compName)
compId = find(strcmpi(model.compNames, compName), 1);
if isempty(compId)
    error('Compartment not found: %s', compName);
end
names = model.metNames(:);
comps = model.metComps(:);
if iscell(names)
    names = cellfun(@char, names, 'UniformOutput', false);
end
idx = find(strcmpi(names, metName) & comps == compId);
if isempty(idx)
    idx = find(strcmpi(names, metName), 1);
end
if isempty(idx)
    error('Metabolite not found: %s (%s)', metName, compName);
end
idx = idx(1);
end
