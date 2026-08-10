%% Production envelopes for tyrosol strains
% mode: 'measured' or 'cumulative'
% addpath('scripts'); tyrosol_envelopes

mode     = 'measured';
plotOnly = false;
medium   = 'Min';
G        = 40;                 % g/L glucose used for yield calc
tyrMW    = 138.164/1000;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root,'scripts'));
home = char(java.lang.System.getProperty('user.home'));
addpath(genpath(fullfile(home,'Documents','CellFactory-ecYeastGEM','code')));
addpath(genpath(fullfile(home,'Documents','ecFactory','code')));
set(0,'DefaultFigureVisible','off');

model = loadEcModel(fullfile(root,'model','ecTyrosol.mat'));
tid   = find(strcmpi(model.rxns,'new_tyrosol_ex'),1);
umap  = loadUsageMap(fullfile(root,'results','candidates_L2.txt'));

dataDir = fullfile(root,'results','data');
figDir  = fullfile(root,'results','figures');
if ~exist(dataDir,'dir'); mkdir(dataDir); end

exp = [];
if strcmpi(mode,'measured')
    xls = findYieldsFile();
    strains = defineByTyrOHStrainTable('measured', xls);
    exp = readFlaskYields(xls, G);
    writeExpCsv(exp, fullfile(dataDir, sprintf('experimental_yields_G%d.csv',G)));
    doNoMito = true;
    tag = 'no_mito';
else
    strains = defineByTyrOHStrainTable('cumulative');
    doNoMito = false;
    tag = medium;
end

m = model;
if doNoMito; m = knockOutMito(m); end
base = setMinGlucose(m, medium);

fprintf('\n=== %s | %s ===\n', mode, tag);
csvFile = fullfile(dataDir, sprintf('envelope_curves_%s_%s.csv', medium, tag));

if plotOnly && exist(csvFile,'file')
    curves = loadCurves(csvFile);
    fprintf('  loaded %s\n', csvFile);
else
    curves = runEnvelopes(base, strains, tid, tyrMW, umap);
    writeCurves(curves, csvFile);
end

if strcmpi(mode,'measured')
    outDir = fullfile(figDir, ['envelopes_' tag]);
    out = fullfile(outDir, sprintf('measured_strain_envelopes_%s_%s', medium, tag));
    cols = strainColors();
else
    outDir = figDir;
    out = fullfile(outDir, sprintf('ByTyrOH_envelopes_%s', medium));
    cols = [];
end
if ~exist(outDir,'dir'); mkdir(outDir); end
plotEnvelopes(curves, exp, out, cols);
fprintf('\nDone.\n');

%% helpers
function m = loadEcModel(f)
raw = load(f); fn = fieldnames(raw); m = raw.(fn{1});
for f = {'lb','ub','c','b','rev'}
    if isfield(m,f{1}); m.(f{1}) = full(double(m.(f{1})(:))); end
end
for f = {'rxns','rxnNames','genes','grRules'}
    if isfield(m,f{1})
        v = m.(f{1}); if isstring(v)||ischar(v); v = cellstr(v); end
        m.(f{1}) = v(:);
    end
end
if isfield(m,'S') && ~issparse(m.S); m.S = sparse(m.S); end
end

function b = setMinGlucose(m, medium)
b = changeMedia_batch(m,'D-glucose exchange (reversible)',medium);
b = setParam(b,'lb',find(strcmpi(b.rxnNames,'growth')),0);
if any(strcmpi(b.rxns,'r_2111'))
    b = setParam(b,'lb','r_2111',0);
    b = setParam(b,'ub','r_2111',1000);
end
end

function m = knockOutMito(m)
% mtDNA genes -> petite
for g = {'Q0045','Q0080','Q0085','Q0105','Q0130','Q0250','Q0275'}
    if any(strcmpi(m.genes,g{1})) || any(strcmpi(m.enzGenes,g{1}))
        m = removeGenes(m,g{1},false,false,false);
    end
end
end

function u = loadUsageMap(f)
u = containers.Map;
if ~exist(f,'file'); return; end
T = readtable(f,'FileType','text','Delimiter','\t');
for i = 1:height(T)
    v = T.maxUsageBio(i);
    if isnan(v) || v <= 0; v = max(T.maxUsage(i),1e-9); end
    u(char(T.genes(i))) = 1.01*v;
end
end

function curves = runEnvelopes(base, strains, tid, mw, umap)
curves = repmat(struct('name','','bio',[],'tyr',[]), numel(strains), 1);
for s = 1:numel(strains)
    mut = applyTyrosolMods(base, strains(s).mods, umap);
    [bio,tyr] = oneEnvelope(mut, tid, 1, mw);
    curves(s).name = strains(s).name;
    curves(s).bio = bio;
    curves(s).tyr = tyr;
    fprintf('  %-12s  %d pts\n', strains(s).name, sum(isfinite(bio)));
end
end

function [bio, tyr] = oneEnvelope(model, target, gur, mw)
gIdx = find(strcmpi(model.rxnNames,'growth'),1);
gluc = find(strcmpi(model.rxnNames,'D-glucose exchange (reversible)'),1);
model = setParam(model,'ub',gluc,1.000001*gur);
model = setParam(model,'lb',gluc,0.999999*gur);
model = setParam(model,'ub',model.rxns{gIdx},1000);
model = setParam(model,'lb',model.rxns{gIdx},0);
model = setParam(model,'ub',model.rxns{target},1000);

sol = solveLP(setParam(model,'obj',model.rxns{gIdx},1));
bio = nan(11,1); tyr = nan(11,1);
if isempty(sol) || ~isfield(sol,'x') || isempty(sol.x) || sol.x(gIdx) < 1e-9
    return;
end

muMax = sol.x(gIdx);
for i = 0:10
    mu = muMax * i / 10;
    m = setParam(model,'obj',model.rxns{target},1);
    m = setParam(m,'lb',m.rxns{gIdx},0.9999*mu);
    sol = solveLP(m,0);
    if isempty(sol) || ~isfield(sol,'x') || isempty(sol.f); continue; end
    q = sol.x(gluc); if abs(q) < 1e-12; continue; end
    bio(i+1) = mu / abs(q * 0.180156);
    tyr(i+1) = sol.x(target) * mw / abs(q * 0.180156);
end
end

function plotEnvelopes(curves, exp, out, cols)
fig = figure('Color','w','Visible','off','Position',[100 100 832 832]);
ax = axes('Parent',fig,'Color','w','LineWidth',1.8,'FontSize',14, ...
    'TickDir','out','TickLength',[0.018 0.025]);
hold(ax,'on'); axis(ax,'square');
maxB = 0; maxT = 0; names = {};

for s = 1:numel(curves)
    ok = isfinite(curves(s).bio) & isfinite(curves(s).tyr);
    if ~any(ok); continue; end
    c = getColor(cols, curves(s).name);
    lw = 2.8; if strcmpi(curves(s).name,'By4743 wt'); lw = 3.4; end
    plot(ax, curves(s).bio, curves(s).tyr, '-', 'LineWidth',lw, 'Color',c, ...
        'DisplayName',curves(s).name);
    names{end+1} = curves(s).name; %#ok<AGROW>
    maxB = max(maxB, max(curves(s).bio(ok)));
    maxT = max(maxT, max(curves(s).tyr(ok)));
    if ~isempty(exp)
        j = find(strcmp({exp.name}, curves(s).name), 1);
        if ~isempty(j)
            scatter(ax, exp(j).bio, exp(j).tyr, 72, 'o', ...
                'MarkerFaceColor',c, 'MarkerEdgeColor','k', 'LineWidth',1.2, ...
                'HandleVisibility','off');
            maxB = max(maxB, max(exp(j).bio));
            maxT = max(maxT, max(exp(j).tyr));
        end
    end
end

legend(ax, names, 'Location','northeast', 'Interpreter','none', 'Box','off', 'FontSize',11);
xlabel(ax,'Biomass yield [gDW/g glucose]');
ylabel(ax,'Tyrosol yield [g/g glucose]');
xlim(ax,[0 max(0.01,1.05*maxB)]);
ylim(ax,[0 max(0.001,1.10*maxT)]);
hold(ax,'off'); box(ax,'on');

savefig(fig,[out '.fig']);
exportgraphics(fig,[out '.png'],'Resolution',600,'BackgroundColor','white');
exportgraphics(fig,[out '.svg'],'ContentType','vector','BackgroundColor','white');
close(fig);
fprintf('  saved %s.{fig,png,svg}\n', out);
end

function exp = readFlaskYields(file, G)
T = readtable(file,'VariableNamingRule','preserve');
sc = findCol(T,{'STRAIN','Strain','strain'});
exp = repmat(struct('name','','bio',[],'tyr',[]), height(T), 1);
for i = 1:height(T)
    od  = [T.(findCol(T,{'OD1'}))(i) T.(findCol(T,{'OD2'}))(i) T.(findCol(T,{'OD3'}))(i)];
    tyr = [T.(findCol(T,{'Tyrosol-1'}))(i) T.(findCol(T,{'Tyrosol-2'}))(i) T.(findCol(T,{'Tyrosol-3'}))(i)];
    exp(i).name = strtrim(char(string(T.(sc)(i))));
    exp(i).bio = (od+0.099)/3.6877/G;
    exp(i).tyr = tyr/(1000*G);
end
end

function writeExpCsv(exp, f)
strain = strings(0); rep = strings(0); bio = []; tyr = [];
for i = 1:numel(exp)
    for k = 1:numel(exp(i).bio)
        strain(end+1,1) = string(exp(i).name); %#ok<AGROW>
        rep(end+1,1) = "R"+k; %#ok<AGROW>
        bio(end+1,1) = exp(i).bio(k); %#ok<AGROW>
        tyr(end+1,1) = exp(i).tyr(k); %#ok<AGROW>
    end
end
writetable(table(strain,rep,bio,tyr,'VariableNames',{'strain','rep','bioYield','tyrYield'}), f);
fprintf('  saved %s\n', f);
end

function writeCurves(curves, f)
strain = strings(0); pt = []; bio = []; tyr = [];
for i = 1:numel(curves)
    n = numel(curves(i).bio);
    strain = [strain; repmat(string(curves(i).name),n,1)]; %#ok<AGROW>
    pt = [pt; (1:n)']; %#ok<AGROW>
    bio = [bio; curves(i).bio(:)]; %#ok<AGROW>
    tyr = [tyr; curves(i).tyr(:)]; %#ok<AGROW>
end
writetable(table(strain,pt,bio,tyr,'VariableNames',{'strain','point','bioYield','tyrYield'}), f);
fprintf('  saved %s\n', f);
end

function curves = loadCurves(f)
T = readtable(f);
curves = repmat(struct('name','','bio',[],'tyr',[]), 0, 1);
for name = unique(T.strain,'stable')'
    rows = strcmp(string(T.strain), string(name));
    curves(end+1).name = char(name); %#ok<AGROW>
    curves(end).bio = T.bioYield(rows);
    curves(end).tyr = T.tyrYield(rows);
end
end

function cols = strainColors()
% gem12; wt black, def gray
gem = [0.1216 0.4667 0.7059; 0.6824 0.7804 0.9098; 0.5961 0.8745 0.5412; ...
       0.8392 0.1529 0.1569; 0.5804 0.4039 0.7412; 0.5490 0.3373 0.2941; ...
       0.8902 0.4667 0.7608; 0.4980 0.4980 0.4980; 0.7373 0.7412 0.1333; ...
       0.0902 0.7451 0.8118; 0.4157 0.2392 0.6039; 0.6941 0.3490 0.1569];
idx = [1 2 3 4 5 6 7 9 10 11];
toh = {'ByTOH1','ByTOH2','ByTOH3','ByTOH4','ByTOH5','ByTOH6','ByTOH7','ByTOH8','ByTOH9','ByTOH10'};
cols = containers.Map('KeyType','char','ValueType','any');
cols('By4743 wt') = [0 0 0];
cols('ByTOHdef') = [0.45 0.45 0.45];
for i = 1:10
    cols(toh{i}) = gem(idx(i),:);
end
end

function rgb = getColor(cols, name)
if isempty(cols); rgb = [0.12 0.20 0.65]; return; end
name = strtrim(char(name));
if isKey(cols,name); rgb = cols(name); else; rgb = [0 0 0]; end
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
for i = 1:numel(names)
    hit = T.Properties.VariableNames(strcmpi(T.Properties.VariableNames, names{i}));
    if ~isempty(hit); c = hit{1}; return; end
end
error('Missing column.');
end
