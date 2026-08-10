function strains = defineByTyrOHStrainTable(mode, yieldsFile)
%DEFINEBYTYROHSTRAINTABLE Genotypes for ByTyrOH / ByTOH envelope simulations.
%
%   strains = defineByTyrOHStrainTable('cumulative')
%       Full cumulative ByTyrOH series (By4743 wt through ByTyrOH def).
%
%   strains = defineByTyrOHStrainTable('measured')
%       Strains listed in the yields workbook only, each with the full
%       cumulative genotype from ByTyrOH_strain_table.md (not single edits).

if nargin < 1 || isempty(mode)
    mode = 'cumulative';
end

switch lower(mode)
    case 'cumulative'
        strains = defineCumulativeByTyrOH();
    case 'measured'
        if nargin < 2 || isempty(yieldsFile)
            yieldsFile = defaultYieldsFile();
        end
        strains = measuredStrainsFromYields(yieldsFile);
    otherwise
        error('Unknown mode: %s (use ''cumulative'' or ''measured'')', mode);
end
end

% -------------------------------------------------------------------------
function strains = defineCumulativeByTyrOH()
k = ecFactoryFactors();

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
strains(9).mods = [strains(8).mods; geneRow('YDR127W', 'OE', k.ARO1oe)];

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

% -------------------------------------------------------------------------
function strains = measuredStrainsFromYields(yieldsFile)
assert(exist(yieldsFile, 'file') == 2, 'Yields file not found: %s', yieldsFile);
opts = detectImportOptions(yieldsFile, 'VariableNamingRule', 'preserve');
T = readtable(yieldsFile, opts);
strainCol = pickColumn(T, {'STRAIN', 'Strain', 'strain'});

cumulative = defineCumulativeByTyrOH();
cumByName = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:numel(cumulative)
    cumByName(cumulative(i).name) = cumulative(i);
end

strains = repmat(struct('name', '', 'mods', {{}}), height(T), 1);
fprintf('Measured strains — full cumulative phenotype (%d):\n', height(T));
for i = 1:height(T)
    yieldName = strtrim(char(string(T.(strainCol)(i))));
    cumKey = cumulativeKeyForYieldName(yieldName);
    if ~isKey(cumByName, cumKey)
        error('No cumulative genotype for %s (key %s)', yieldName, cumKey);
    end
    entry = cumByName(cumKey);
    strains(i).name = yieldName;
    strains(i).mods = entry.mods;
    fprintf('  %-12s  <- %-12s  (%d edit(s))\n', yieldName, cumKey, size(entry.mods, 1));
end
end

function key = cumulativeKeyForYieldName(yieldName)
% Map flask labels (ByTOH1…) to cumulative table (ByTyrOH 1…; step 4 omitted).
n = lower(strtrim(char(yieldName)));
if strcmp(n, 'by4743 wt')
    key = 'By4743 wt';
    return;
end
if strcmp(n, 'bytohdef')
    key = 'ByTyrOH def';
    return;
end
tok = regexp(n, '^bytoh(\d+)$', 'tokens', 'once');
if isempty(tok)
    error('Unknown strain name in yields file: %s', yieldName);
end
idx = str2double(tok{1});
cumStep = [1 2 3 5 6 7 8 9 10 11];
if idx < 1 || idx > numel(cumStep)
    error('ByTOH index out of range: %s', yieldName);
end
key = sprintf('ByTyrOH %d', cumStep(idx));
end

% -------------------------------------------------------------------------
function mods = modificationToMods(modStr, k)
mods = cell(0, 3);
mod = lower(strtrim(char(modStr)));
if isempty(mod) || strcmp(mod, 'none') || strcmp(mod, 'nan')
    return;
end

if contains(mod, 'aro2')
    mods = [mods; geneRow('YGL148W', 'OE', k.ARO2)]; %#ok<AGROW>
end
if contains(mod, 'aro7')
    mods = [mods; geneRow('YPR060C', 'OE', k.ARO7)]; %#ok<AGROW>
end
if contains(mod, 'aro4')
    mods = [mods; geneRow('YBR249C', 'OE', k.ARO4)]; %#ok<AGROW>
end
if contains(mod, 'aro10')
    mods = [mods; geneRow('YDR380W', 'OE', k.ARO10)]; %#ok<AGROW>
end
if contains(mod, 'lpp1')
    mods = [mods; geneRow('YDR503C', 'KD', k.LPP1)]; %#ok<AGROW>
end
if contains(mod, 'ald5')
    mods = [mods; geneRow('YER073W', 'KO', 0)]; %#ok<AGROW>
end
if contains(mod, 'mae1')
    mods = [mods; geneRow('YKL029C', 'KO', 0)]; %#ok<AGROW>
end
if contains(mod, 'zwf1')
    mods = [mods; geneRow('YNL241C', 'KD', k.ZWF1)]; %#ok<AGROW>
end
if contains(mod, 'pha2')
    mods = [mods; geneRow('YNL316C', 'KD', k.PHA2)]; %#ok<AGROW>
end
if contains(mod, 'adh6')
    mods = [mods; geneRow('YMR318C', 'OE', k.ADH6)]; %#ok<AGROW>
end
if contains(mod, 'adh7')
    mods = [mods; geneRow('YCR105W', 'OE', k.ADH7)]; %#ok<AGROW>
end
if contains(mod, 'pdh1') || contains(mod, 'mtpdh')
    mods = [mods; geneRow('YER178W', 'OE', k.PDH1)]; %#ok<AGROW>
end
if contains(mod, 'aro1')
    if contains(mod, 'kd')
        mods = [mods; geneRow('YDR127W', 'KD', k.ARO1kd)]; %#ok<AGROW>
    else
        mods = [mods; geneRow('YDR127W', 'OE', k.ARO1oe)]; %#ok<AGROW>
    end
end

if isempty(mods)
    error('Could not parse modification string: %s', modStr);
end
end

% -------------------------------------------------------------------------
function k = ecFactoryFactors()
k = struct( ...
    'ARO2', 13.9306, 'ARO7', 15.8562, 'ARO4', 13.9306, 'ARO10', 1000, ...
    'ARO1oe', 13.9306, 'ADH6', 13.9306, 'ADH7', 1000, 'PDH1', 5.2514, ...
    'LPP1', 0.2079, 'PHA2', 0.21, 'ZWF1', 0.21, 'ARO1kd', 0.21);
end

% -------------------------------------------------------------------------
function f = defaultYieldsFile()
homeDir = char(java.lang.System.getProperty('user.home'));
candidates = { ...
    fullfile(homeDir, 'Library', 'CloudStorage', 'OneDrive-Chalmers', ...
        'Documents', 'tyrosol_ecYeasy', '6. Cepas tirosol - rendimientos en matraz.xlsx'), ...
    fullfile(homeDir, 'Documents', 'tyrosol_ecYeasy', ...
        '6. Cepas tirosol - rendimientos en matraz.xlsx')};
for i = 1:numel(candidates)
    if exist(candidates{i}, 'file')
        f = candidates{i};
        return;
    end
end
error('Yields workbook not found. Pass yieldsFile explicitly.');
end

function col = pickColumn(T, candidates)
vars = T.Properties.VariableNames;
for i = 1:numel(candidates)
    hit = vars(strcmpi(vars, candidates{i}));
    if ~isempty(hit)
        col = hit{1};
        return;
    end
end
error('Expected one of these columns: %s', strjoin(candidates, ', '));
end

function s = modStrForRow(v)
if iscell(v); s = char(v{1}); else; s = char(string(v)); end
end

function row = geneRow(gene, action, factor)
row = {gene, action, factor};
end
