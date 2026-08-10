function build_ecTyrosol_model_raven(outFile, baseFile)
%BUILD_ECTYROSOL_MODEL_RAVEN Build ecTyrosol.mat with RAVEN Toolbox.
%
%   Extends ecYeastGEM_batch.mat for tyrosol via the native Ehrlich pathway
%   using standard RAVEN functions (addMets, addRxns, setParam). No direct
%   hand-editing of the saved .mat file.
%
%   Assumptions A1–A6: docs/METHODS.md
%   Enzyme arms: kcat = 1000 1/s on ARO10 and ADH7 (coef = 1/(kcat*3600)).
%
%   Usage:
%     cd model
%     build_ecTyrosol_model_raven
%
%   Then strain design:
%     addpath('../scripts'); run_ecTyrosol_native
%
%   Optional environment variable: ECYEASTGEM_BATCH (path to base .mat)
%
%   Requires RAVEN Toolbox on the MATLAB path.

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
assert(~isempty(which('addRxns')), 'RAVEN addRxns not on path.');
assert(~isempty(which('addMets')), 'RAVEN addMets not on path.');

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

model.c(:) = 0;

metsToAdd.mets = {'s_4hpaa_c', 's_tyrosol_c', 's_tyrosol_e'};
metsToAdd.metNames = {'4-hydroxyphenylacetaldehyde', 'tyrosol', 'tyrosol'};
metsToAdd.compartments = {'c', 'c', 'e'};
if isfield(model, 'metFormulas')
    metsToAdd.metFormulas = {'C8H8O2', 'C8H10O2', 'C8H10O2'};
end
model = addMets(model, metsToAdd);

mH = metId(model, 's_0794');
mCO2 = metId(model, 's_0456');
mNADPH = metId(model, 's_1212');
mNADP = metId(model, 's_1207');
mHpp = metId(model, findMetByNameComp(model, '3-(4-hydroxyphenyl)pyruvate', 'cytoplasm'));
mProtAro10 = metId(model, 'prot_Q06408');
mProtAdh7 = metId(model, 'prot_P25377');

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

if exist('setParam', 'file')
    model = setParam(model, 'obj', 'new_tyrosol_ex', 1);
else
    model.c(:) = 0;
    model.c(strcmp(model.rxns, 'new_tyrosol_ex')) = 1;
end

if isfield(model, 'description')
    model.description = [ ...
        'ecTyrosol: ecYeastGEM_batch + Ehrlich tyrosol pathway (ARO10, ADH7). ', ...
        'Built with RAVEN addMets/addRxns.'];
end

[outDir, ~, ~] = fileparts(outFile);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end

out = struct();
out.(varName) = model;
save(outFile, '-struct', 'out', '-v7');
fprintf('Saved %s (%d rxns, %d mets)\n', outFile, numel(model.rxns), numel(model.mets));

nativeAlias = fullfile(outDir, 'ecTyrosol_native.mat');
if ~strcmp(nativeAlias, outFile)
    copyfile(outFile, nativeAlias);
    fprintf('Copied alias: %s\n', nativeAlias);
end

fprintf('Product objective: new_tyrosol_ex\n');
fprintf('Next: addpath(''../scripts''); run_ecTyrosol_native\n');

end

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
