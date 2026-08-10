function strains = defineByTyrOHStrainTable(mode, xls)
% Cumulative genotypes for ByTyrOH / ByTOH strains.
% mode = 'cumulative' | 'measured'

if nargin < 1 || isempty(mode); mode = 'cumulative'; end

switch lower(mode)
    case 'cumulative'
        strains = byTyrOH();
    case 'measured'
        if nargin < 2 || isempty(xls); xls = findYieldsFile(); end
        strains = fromYields(xls);
    otherwise
        error('mode must be cumulative or measured');
end
end

function s = byTyrOH()
% step 4 skipped (lab numbering)
k = factors();

s(1).name = 'By4743 wt';  s(1).mods = cell(0,3);
s(2).name = 'ByTyrOH 1';  s(2).mods = row('YGL148W','OE',k.ARO2);
s(3).name = 'ByTyrOH 2';  s(3).mods = [s(2).mods; row('YPR060C','OE',k.ARO7)];
s(4).name = 'ByTyrOH 3';  s(4).mods = [s(3).mods; row('YBR249C','OE',k.ARO4); row('YDR380W','OE',k.ARO10)];
s(5).name = 'ByTyrOH 5';  s(5).mods = [s(4).mods; row('YDR503C','KD',k.LPP1)];
s(6).name = 'ByTyrOH 6';  s(6).mods = [s(5).mods; row('YER073W','KO',0)];
s(7).name = 'ByTyrOH 7';  s(7).mods = [s(6).mods; row('YKL029C','KO',0)];
s(8).name = 'ByTyrOH 8';  s(8).mods = [s(7).mods; row('YNL241C','KD',k.ZWF1)];
s(9).name = 'ByTyrOH 9';  s(9).mods = [s(8).mods; row('YDR127W','OE',k.ARO1)];
s(10).name = 'ByTyrOH 10'; s(10).mods = [s(9).mods; row('YNL316C','KD',k.PHA2)];
s(11).name = 'ByTyrOH 11'; s(11).mods = [s(10).mods; row('YMR318C','OE',k.ADH6); row('YCR105W','OE',k.ADH7)];
s(12).name = 'ByTyrOH 12'; s(12).mods = [s(10).mods; row('YER178W','OE',k.PDH1)];
s(13).name = 'ByTyrOH def'; s(13).mods = [s(11).mods; row('YER178W','OE',k.PDH1)];
end

function s = fromYields(xls)
assert(exist(xls,'file')==2, 'File not found: %s', xls);
opts = detectImportOptions(xls,'VariableNamingRule','preserve');
T = readtable(xls, opts);
col = findCol(T, {'STRAIN','Strain','strain'});

all = byTyrOH();
map = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(all); map(all(i).name) = all(i); end

s = repmat(struct('name','','mods',{{}}), height(T), 1);
fprintf('Measured strains (%d):\n', height(T));
for i = 1:height(T)
    name = strtrim(char(string(T.(col)(i))));
    key  = mapName(name);
    if ~isKey(map,key); error('No genotype for %s', name); end
    s(i).name = name;
    s(i).mods = map(key).mods;
    fprintf('  %-12s <- %s (%d)\n', name, key, size(s(i).mods,1));
end
end

function key = mapName(name)
% ByTOH1..11 -> ByTyrOH 1,2,3,5..11
n = lower(strtrim(char(name)));
if strcmp(n,'by4743 wt'); key = 'By4743 wt'; return; end
if strcmp(n,'bytohdef');  key = 'ByTyrOH def'; return; end
tok = regexp(n,'^bytoh(\d+)$','tokens','once');
if isempty(tok); error('Unknown strain: %s', name); end
idx = str2double(tok{1});
steps = [1 2 3 5 6 7 8 9 10 11];
if idx < 1 || idx > numel(steps); error('Bad ByTOH index: %s', name); end
key = sprintf('ByTyrOH %d', steps(idx));
end

function k = factors()
k = struct( ...
    'ARO2',13.9306,'ARO7',15.8562,'ARO4',13.9306,'ARO10',1000, ...
    'ARO1',13.9306,'ADH6',13.9306,'ADH7',1000,'PDH1',5.2514, ...
    'LPP1',0.2079,'PHA2',0.21,'ZWF1',0.21);
end

function f = findYieldsFile()
h = char(java.lang.System.getProperty('user.home'));
c = {fullfile(h,'Library','CloudStorage','OneDrive-Chalmers','Documents','tyrosol_ecYeasy','6. Cepas tirosol - rendimientos en matraz.xlsx'), ...
     fullfile(h,'Documents','tyrosol_ecYeasy','6. Cepas tirosol - rendimientos en matraz.xlsx')};
for i = 1:2
    if exist(c{i},'file'); f = c{i}; return; end
end
error('Yields Excel not found.');
end

function c = findCol(T, names)
v = T.Properties.VariableNames;
for i = 1:numel(names)
    hit = v(strcmpi(v,names{i}));
    if ~isempty(hit); c = hit{1}; return; end
end
error('Missing column.');
end

function r = row(gene, act, fac)
r = {gene, act, fac};
end
