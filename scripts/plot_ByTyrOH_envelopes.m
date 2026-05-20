function plot_ByTyrOH_envelopes(modelFile, mediumName, outDir)
%PLOT_BYTYROH_ENVELOPES Tyrosol production envelopes for ByTyrOH lab strains.
%
%   Uses the same ecModel as the Iván-matching ecFactory run:
%     model/ecTyrosol.mat  (native Ehrlich pathway, ARO10/ADH7 kcat = 1000)
%     product exchange: new_tyrosol_ex
%     medium: Min by default (matches run_tyrosol_ecFactory.m; change if needed)
%
%   One biomass–tyrosol envelope per cumulative mutant (getYieldPlot +
%   getMutantModel), all curves in one figure.
%
%   Usage:
%     addpath('scripts')
%     plot_ByTyrOH_envelopes
%     plot_ByTyrOH_envelopes([], 'YEP')   % optional medium override
%
%   The reference single-curve figure (e.g. ecTyrosol_fma on YEP) is only an
%   example of the plot type; this script uses ecTyrosol.mat for Iván targets.
%
%   Output: results/figures/ByTyrOH_envelopes_<medium>.png

pkgRoot = fileparts(fileparts(mfilename('fullpath')));

if nargin < 1 || isempty(modelFile)
    modelFile = fullfile(pkgRoot, 'model', 'ecTyrosol.mat');
end
if nargin < 2 || isempty(mediumName)
    mediumName = 'Min';
end
if nargin < 3 || isempty(outDir)
    outDir = fullfile(pkgRoot, 'results', 'figures');
end

homeDir = char(java.lang.System.getProperty('user.home'));
addpath(genpath(fullfile(homeDir, 'Documents', 'CellFactory-ecYeastGEM', 'code')));
addpath(genpath(fullfile(homeDir, 'Documents', 'ecFactory', 'code')));

assert(exist(modelFile, 'file') == 2, ...
    'ecTyrosol.mat not found at %s. Run build_ecTyrosol_model_raven.m first.', modelFile);
if ~exist(outDir, 'dir'); mkdir(outDir); end

raw = load(modelFile);
fn = fieldnames(raw);
if numel(fn) == 1
    ecModel = prepareEcModel(raw.(fn{1}));
else
    load(modelFile, 'model');
    ecModel = prepareEcModel(model);
end

targetRxn = 'new_tyrosol_ex';
targetIdx = find(strcmpi(ecModel.rxns, targetRxn), 1);
assert(~isempty(targetIdx), ...
    'Target %s not in %s. This script requires the native-Ehrlich tyrosol model.', ...
    targetRxn, modelFile);

carbonSource = 'D-glucose exchange (reversible)';
tyrosolMW = 138.164 / 1000;
GUR = 1;

baseModel = changeMedia_batch(ecModel, carbonSource, mediumName);
baseModel = setParam(baseModel, 'lb', find(strcmpi(baseModel.rxnNames, 'growth')), 0);
if any(strcmpi(baseModel.rxns, 'r_2111'))
    baseModel = setParam(baseModel, 'lb', 'r_2111', 0);
    baseModel = setParam(baseModel, 'ub', 'r_2111', 1000);
end
usageMap = loadCandidateEnzUsage(fullfile(pkgRoot, 'results', 'candidates_L2.txt'));
strains = defineByTyrOHStrains();
nStrains = numel(strains);

fprintf('=== ByTyrOH production envelopes ===\n');
fprintf('Model : %s (native Ehrlich, Iván-matching build)\n', modelFile);
fprintf('Medium: %s | Target: %s | GUR: %g\n', mediumName, targetRxn, GUR);

fig = figure('Color', 'w', 'Position', [100 100 960 640]);
hold on;
legends = cell(nStrains, 1);
maxBio = 0;
maxTyr = 0;

for s = 1:nStrains
    strain = strains(s);
    mutModel = applyStrainMods(baseModel, strain.mods, usageMap);
    [bioYield, tyrYield] = getYieldPlotTyrosol(mutModel, targetIdx, GUR, tyrosolMW);
    if all(isnan(bioYield))
        warning('Strain %s: no envelope points; skipping.', strain.name);
        legends{s} = '';
        continue;
    end
    maxBio = max(maxBio, max(bioYield(~isnan(bioYield))));
    maxTyr = max(maxTyr, max(tyrYield(~isnan(tyrYield))));

    if s == 1
        plot(bioYield, tyrYield, '-', 'LineWidth', 2.5, 'Color', [0.12 0.20 0.65]);
    else
        plot(bioYield, tyrYield, '-', 'LineWidth', 1.8);
    end
    legends{s} = strain.name;
    fprintf('  %-12s  max tyr = %.4e g/g  max bio = %.4f  (%d/%d points)\n', ...
        strain.name, max(tyrYield(~isnan(tyrYield))), max(bioYield(~isnan(bioYield))), ...
        sum(~isnan(bioYield)), numel(bioYield));
end

hold off;
box on;
xlabel('Biomass yield [gDW/g glucose]');
ylabel('Tyrosol yield [g/g glucose]');
title(sprintf('ByTyrOH envelopes (ecTyrosol, %s)', mediumName));
legends = legends(~cellfun(@isempty, legends));
legend(legends, 'Location', 'northeast', 'Interpreter', 'none');
set(gca, 'FontSize', 12);
xlim([0, max(0.01, 1.05 * maxBio)]);
ylim([0, max(0.001, 1.10 * maxTyr)]);

tag = mediumName;
pngFile = fullfile(outDir, sprintf('ByTyrOH_envelopes_%s.png', tag));
figFile = fullfile(outDir, sprintf('ByTyrOH_envelopes_%s.fig', tag));
saveas(fig, pngFile);
saveas(fig, figFile);
fprintf('Saved:\n  %s\n  %s\n', pngFile, figFile);

end

% -------------------------------------------------------------------------
function mutModel = applyStrainMods(baseModel, mods, usageMap)
% Cumulative ecFactory edits (k_score factors + maxUsageBio from candidates).
mutModel = baseModel;
for m = 1:size(mods, 1)
    fac = mods{m, 3};
    if iscell(fac); fac = fac{1}; end
    gene = char(mods{m, 1});
    action = char(mods{m, 2});
    if shouldSkipNativeMod(mutModel, gene, action)
        continue;
    end
    edit = {gene, action, double(fac)};
    baseUsage = lookupEnzUsage(usageMap, edit{1}, mutModel);
    mutModel = applyEcFactoryMod(mutModel, edit, baseUsage);
end
end

function skip = shouldSkipNativeMod(model, gene, action)
% Native ecTyrosol already wires ARO10/ADH7 at kcat=1000 on new pathway rxns.
skip = false;
if ~strcmpi(action, 'OE')
    return;
end
if strcmpi(gene, 'YDR380W') && any(strcmp(model.rxns, 'new_aro10_HPP'))
    skip = true;
end
if strcmpi(gene, 'YCR105W') && any(strcmp(model.rxns, 'new_adh7_tyrosol'))
    skip = true;
end
end

function usageMap = loadCandidateEnzUsage(candFile)
usageMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
if ~exist(candFile, 'file')
    return;
end
T = readtable(candFile, 'FileType', 'text', 'Delimiter', '\t');
for i = 1:height(T)
    g = char(T.genes(i));
    u = T.maxUsageBio(i);
    if isnan(u) || u <= 0
        u = max(T.maxUsage(i), 1e-9);
    end
    usageMap(g) = 1.01 * u;
end
end

function baseUsage = lookupEnzUsage(usageMap, gene, model)
if nargin > 1 && ~isempty(usageMap) && isKey(usageMap, gene)
    baseUsage = usageMap(gene);
    return;
end
baseUsage = enzymeBaseUsage(model, gene);
end

function mutantModel = applyEcFactoryMod(model, modifications, base_usage)
% ecFactory getMutantModel logic without the post-hoc viability error (line 117).
mutantModel = model;
genes2mod = modifications(:, 1);
actions = modifications(:, 2);
expF = modifications(:, 3);

for i = 1:numel(genes2mod)
    gene = genes2mod{i};
    action = actions{i};
    expFactor = expF{i};
    gene2modIndex = findEnzGeneIndex(mutantModel, gene);
    if strcmpi(action, 'KO')
        mutantModel = removeGenes(mutantModel, gene, false, false, false);
    end
    if ~isempty(gene2modIndex)
        enzyme = mutantModel.enzymes{gene2modIndex(1)};
        if strcmpi(action, 'KD') || strcmpi(action, 'OE')
            enzRxn = find(contains(mutantModel.rxnNames, enzyme));
            for j = 1:numel(enzRxn)
                r = enzRxn(j);
                if contains(mutantModel.rxnNames{r}, 'exchange') || mutantModel.ub(r) < 100
                    mutantModel.ub(r) = base_usage * expFactor;
                    mutantModel.lb(r) = 0;
                elseif ~isempty(base_usage)
                    if expFactor > 1
                        mutantModel.lb(r) = min(1000, base_usage * expFactor);
                        mutantModel.ub(r) = 1000;
                    else
                        mutantModel.ub(r) = base_usage * expFactor;
                        mutantModel.lb(r) = 0;
                    end
                end
                if mutantModel.ub(r) <= mutantModel.lb(r)
                    mutantModel.lb(r) = 0.99 * mutantModel.ub(r);
                end
            end
        end
    elseif strcmpi(action, 'KD')
        geneRxns = find(contains(mutantModel.grRules, gene));
        if ~isempty(geneRxns)
            mutantModel.S(:, geneRxns) = mutantModel.S(:, geneRxns) * expFactor;
        end
    end
end
end

function idx = findEnzGeneIndex(model, gene)
idx = find(strcmpi(model.enzGenes, gene), 1);
if ~isempty(idx)
    return;
end
for i = 1:numel(model.enzGenes)
    entry = model.enzGenes{i};
    if iscell(entry)
        entry = entry{1};
    end
    if ischar(entry) || isstring(entry)
        if strcmpi(char(entry), gene)
            idx = i;
            return;
        end
    end
end
idx = [];
end

% -------------------------------------------------------------------------
function ecModel = prepareEcModel(ecModel)
vectorFields = {'lb','ub','c','b','rev','metCharges','metComps','rxnConfidenceScores'};
for i = 1:numel(vectorFields)
    field = vectorFields{i};
    if isfield(ecModel, field)
        ecModel.(field) = full(double(ecModel.(field)(:)));
    end
end
cellFields = {'rxns','rxnNames','mets','metNames','genes','grRules','rules', ...
              'eccodes','rxnReferences','subSystems','metFormulas'};
for i = 1:numel(cellFields)
    field = cellFields{i};
    if isfield(ecModel, field)
        value = ecModel.(field);
        if ischar(value); value = cellstr(value);
        elseif isstring(value); value = cellstr(value);
        end
        ecModel.(field) = value(:);
    end
end
if isfield(ecModel, 'S') && ~issparse(ecModel.S)
    ecModel.S = sparse(ecModel.S);
end
end

% -------------------------------------------------------------------------
function strains = defineByTyrOHStrains()
% Cumulative ByTyrOH table. OE/KD factors = ecFactory k_scores on ecTyrosol (Iván L3).
k = struct( ...
    'ARO2', 13.9306, 'ARO7', 15.8562, 'ARO4', 13.9306, 'ARO10', 1000, ...
    'ARO1', 13.9306, 'ADH6', 13.9306, 'ADH7', 1000, 'PDH1', 5.2514, ...
    'LPP1', 0.2079, 'PHA2', 0.21, 'ZWF1', 0.21);

strains(1).name = 'By4743 wt';
strains(1).mods = cell(0, 3);

strains(2).name = 'ByTyrOH 1';
strains(2).mods = geneRow('YGL148W', 'OE', k.ARO2);

strains(3).name = 'ByTyrOH 2';
strains(3).mods = [strains(2).mods; geneRow('YPR060C', 'OE', k.ARO7)];

strains(4).name = 'ByTyrOH 3';
strains(4).mods = [strains(3).mods; ...
    geneRow('YBR249C', 'OE', k.ARO4); ...
    geneRow('YDR380W', 'OE', k.ARO10)];

strains(5).name = 'ByTyrOH 5';
strains(5).mods = [strains(4).mods; geneRow('YDR503C', 'KD', k.LPP1)];

strains(6).name = 'ByTyrOH 6';
strains(6).mods = [strains(5).mods; geneRow('YER073W', 'KO', 0)];

strains(7).name = 'ByTyrOH 7';
strains(7).mods = [strains(6).mods; geneRow('YKL029C', 'KO', 0)];

strains(8).name = 'ByTyrOH 8';
strains(8).mods = [strains(7).mods; geneRow('YNL241C', 'KD', k.ZWF1)];

strains(9).name = 'ByTyrOH 9';
strains(9).mods = [strains(8).mods; geneRow('YDR127W', 'OE', k.ARO1)];

strains(10).name = 'ByTyrOH 10';
strains(10).mods = [strains(9).mods; geneRow('YNL316C', 'KD', k.PHA2)];

strains(11).name = 'ByTyrOH 11';
strains(11).mods = [strains(10).mods; ...
    geneRow('YMR318C', 'OE', k.ADH6); ...
    geneRow('YCR105W', 'OE', k.ADH7)];

strains(12).name = 'ByTyrOH 12';
strains(12).mods = [strains(10).mods; geneRow('YER178W', 'OE', k.PDH1)];

strains(13).name = 'ByTyrOH def';
strains(13).mods = [strains(11).mods; geneRow('YER178W', 'OE', k.PDH1)];
end

function row = geneRow(gene, action, factor)
row = {gene, action, factor};
end

function miu = growthFluxMax(model)
gIdx = find(strcmpi(model.rxnNames, 'growth'), 1);
tmp = setParam(model, 'obj', model.rxns{gIdx}, 1);
sol = solveLP(tmp);
if isempty(sol) || ~isfield(sol, 'x') || isempty(sol.x) || sol.x(gIdx) < 1e-9
    miu = 0.48;
else
    miu = sol.x(gIdx);
end
end

function [BioYield, yield] = getYieldPlotTyrosol(model, target, GUR, MW)
% CellFactory getYieldPlot: each strain uses its own max growth (not WT).
% Failed scan steps are NaN so plot() does not draw spurious connecting segments.
ecM = any(contains(model.rxnNames, 'prot_pool'));
growthIndex = find(strcmpi(model.rxnNames, 'growth'), 1);
if GUR > 1
    lower_uptk = 0;
else
    lower_uptk = GUR;
end
if ~ecM
    glucIndex = find(strcmpi(model.rxnNames, 'D-glucose exchange'), 1);
    model = setParam(model, 'lb', glucIndex, -1.000001 * GUR);
    model = setParam(model, 'ub', glucIndex, -0.999999 * lower_uptk);
else
    glucIndex = find(strcmpi(model.rxnNames, 'D-glucose exchange (reversible)'), 1);
    model = setParam(model, 'ub', glucIndex, 1.000001 * GUR);
    model = setParam(model, 'lb', glucIndex, 0.999999 * lower_uptk);
end
model = setParam(model, 'ub', model.rxns{growthIndex}, 1000);
model = setParam(model, 'lb', model.rxns{growthIndex}, 0);
model = setParam(model, 'ub', model.rxns{target}, 1000);

MiuMax = growthFluxMax(model);
if MiuMax < 1e-9
    BioYield = nan(11, 1);
    yield = nan(11, 1);
    return;
end

iterations = 10;
BioYield = nan(iterations + 1, 1);
yield = nan(iterations + 1, 1);
for i = 1:iterations + 1
    tempModel = setParam(model, 'obj', model.rxns{target}, 1);
    Drate = MiuMax * (i - 1) / iterations;
    tempModel = setParam(tempModel, 'lb', tempModel.rxns{growthIndex}, 0.9999 * Drate);
    tempModel = setParam(tempModel, 'ub', tempModel.rxns{growthIndex}, 1000);
    solution = solveLP(tempModel, 1);
    if isempty(solution) || ~isfield(solution, 'x') || isempty(solution.f)
        continue;
    end
    GURsim = solution.x(glucIndex);
    if abs(GURsim) < 1e-12
        continue;
    end
    production = solution.x(target);
    BioYield(i) = Drate / abs(GURsim * 0.180156);
    yield(i) = production * MW / abs(GURsim * 0.180156);
end
end

function usage = enzymeBaseUsage(model, gene)
idx = find(strcmpi(model.enzGenes, gene), 1);
if isempty(idx); usage = 1e-9; return; end
enzyme = model.enzymes{idx};
enzRxn = find(contains(model.rxnNames, enzyme), 1);
if isempty(enzRxn); usage = 1e-9; return; end
growthIdx = find(strcmpi(model.rxnNames, 'growth'), 1);
temp = setParam(model, 'obj', model.rxns{growthIdx}, 1);
sol = solveLP(temp);
if isempty(sol) || ~isfield(sol, 'x') || isempty(sol.x)
    usage = 1e-9;
else
    usage = max(sol.x(enzRxn), 1e-9);
end
end
