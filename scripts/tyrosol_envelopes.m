%% Tyrosol production envelopes — one script
%  mode = 'measured'  -> 12 flask strains + lab yield dots + with/no mito
%  mode = 'cumulative' -> full ByTyrOH cumulative series, simulation plot only
%
%  addpath('scripts'); tyrosol_envelopes

mode = 'measured';   % 'measured' | 'cumulative'
plotOnly = false;    % set true to replot from CSV only (after full run)
medium = 'Min';
G = 40;              % g/L glucose consumed (experimental yields only)
tyrMW = 138.164 / 1000;

pkg = fileparts(fileparts(mfilename('fullpath')));
addpath(pkg); addpath(fullfile(pkg,'scripts'));
home = char(java.lang.System.getProperty('user.home'));
addpath(genpath(fullfile(home,'Documents','CellFactory-ecYeastGEM','code')));
addpath(genpath(fullfile(home,'Documents','ecFactory','code')));
set(0,'DefaultFigureVisible','off');

m = loadModel(fullfile(pkg,'model','ecTyrosol.mat'));
tid = find(strcmpi(m.rxns,'new_tyrosol_ex'),1);
umap = usageMap(fullfile(pkg,'results','candidates_L2.txt'));
dataDir = fullfile(pkg,'results','data');
figDir = fullfile(pkg,'results','figures');
if ~exist(dataDir,'dir'); mkdir(dataDir); end

exp = [];
if strcmpi(mode,'measured')
    yieldsFile = flaskFile();
    strains = defineByTyrOHStrainTable('measured', yieldsFile);
    exp = readFlaskYields(yieldsFile, G);
    writeExpCsv(exp, fullfile(dataDir, sprintf('experimental_yields_G%d.csv',G)));
    scenarios = true;            % no-mito (rho0) background only
    tags = {'no_mito'};
else
    strains = defineByTyrOHStrainTable('cumulative');
    scenarios = false;
    tags = {medium};
end

for k = 1:numel(scenarios)
    mk = m;
    if scenarios(k); mk = noMito(mk); end
    base = baseModel(mk, medium);
    tag = tags{k};
    fprintf('\n=== %s | %s ===\n', mode, tag);

    envCsv = fullfile(dataDir, sprintf('envelope_curves_%s_%s.csv', medium, tag));
    if plotOnly && exist(envCsv, 'file')
        curves = loadCurveCsv(envCsv);
        fprintf('  (plot only — loaded %s)\n', envCsv);
    else
        curves = simulateEnvelopes(base, strains, tid, tyrMW, umap);
        writeCurveCsv(curves, envCsv);
    end
    if strcmpi(mode,'measured')
        outDir = fullfile(figDir, ['envelopes_' tag]);
    else
        outDir = figDir;
    end
    if ~exist(outDir,'dir'); mkdir(outDir); end
    if strcmpi(mode,'measured')
        out = fullfile(outDir, sprintf('measured_strain_envelopes_%s_%s',medium,tag));
        plotEnvelopes(curves, exp, out, strainColors());
    else
        out = fullfile(outDir, sprintf('ByTyrOH_envelopes_%s',medium));
        plotEnvelopes(curves, [], out, []);
    end
end
fprintf('\nDone.\n');

%% --- local functions ---
function m = loadModel(f)
raw = load(f); fn = fieldnames(raw); m = raw.(fn{1});
vecFields = {'lb','ub','c','b','rev'};
for k = 1:numel(vecFields)
    fld = vecFields{k};
    if isfield(m, fld); m.(fld) = full(double(m.(fld)(:))); end
end
cellFields = {'rxns','rxnNames','genes','grRules'};
for k = 1:numel(cellFields)
    fld = cellFields{k};
    if isfield(m, fld)
        v = m.(fld);
        if isstring(v) || ischar(v); v = cellstr(v); end
        m.(fld) = v(:);
    end
end
if isfield(m,'S') && ~issparse(m.S); m.S = sparse(m.S); end
end

function b = baseModel(m, medium)
b = changeMedia_batch(m,'D-glucose exchange (reversible)',medium);
b = setParam(b,'lb',find(strcmpi(b.rxnNames,'growth')),0);
if any(strcmpi(b.rxns,'r_2111')); b=setParam(b,'lb','r_2111',0); b=setParam(b,'ub','r_2111',1000); end
end

function m = noMito(m)
for g = {'Q0045','Q0080','Q0085','Q0105','Q0130','Q0250','Q0275'}
    if any(strcmpi(m.genes,g{1}))||any(strcmpi(m.enzGenes,g{1})); m=removeGenes(m,g{1},false,false,false); end
end
end

function u = usageMap(f)
u = containers.Map; if ~exist(f,'file'); return; end
T = readtable(f,'FileType','text','Delimiter','\t');
for i=1:height(T); v=T.maxUsageBio(i); if isnan(v)||v<=0; v=max(T.maxUsage(i),1e-9); end; u(char(T.genes(i)))=1.01*v; end
end

function curves = simulateEnvelopes(base, strains, tid, mw, umap)
curves = repmat(struct('name','','bio',[],'tyr',[]),numel(strains),1);
for s=1:numel(strains)
    [bio,tyr] = envelopeCurve(applyTyrosolMods(base,strains(s).mods,umap), tid, 1, mw);
    curves(s).name=strains(s).name; curves(s).bio=bio; curves(s).tyr=tyr;
    fprintf('  %-12s  %d points\n', strains(s).name, sum(isfinite(bio)));
end
end

function [bio, tyr] = envelopeCurve(model, target, GUR, MW)
gIdx = find(strcmpi(model.rxnNames,'growth'),1);
gluc = find(strcmpi(model.rxnNames,'D-glucose exchange (reversible)'),1);
model = setParam(model,'ub',gluc,1.000001*GUR);
model = setParam(model,'lb',gluc,0.999999*GUR);
model = setParam(model,'ub',model.rxns{gIdx},1000);
model = setParam(model,'lb',model.rxns{gIdx},0);
model = setParam(model,'ub',model.rxns{target},1000);
sol = solveLP(setParam(model,'obj',model.rxns{gIdx},1));
bio = nan(11,1); tyr = nan(11,1);
if isempty(sol) || ~isfield(sol,'x') || isempty(sol.x) || sol.x(gIdx) < 1e-9; return; end
miu = sol.x(gIdx);
for i = 0:10
    D = miu * i / 10;
    m = setParam(model,'obj',model.rxns{target},1);
    m = setParam(m,'lb',m.rxns{gIdx},0.9999*D);
    sol = solveLP(m,0);
    if isempty(sol) || ~isfield(sol,'x') || isempty(sol.f); continue; end
    g = sol.x(gluc); if abs(g) < 1e-12; continue; end
    bio(i+1) = D / abs(g * 0.180156);
    tyr(i+1) = sol.x(target) * MW / abs(g * 0.180156);
end
end

function plotEnvelopes(curves, exp, out, cmap)
figW = round(1040 * 0.8);
fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 figW figW]);
ax = axes('Parent', fig, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
    'LineWidth', 1.8, 'FontSize', 14, 'TickDir', 'out', 'TickLength', [0.018 0.025]);
hold(ax, 'on'); axis(ax, 'square');
maxB = 0; maxT = 0; leg = {};
for s = 1:numel(curves)
    if ~any(isfinite(curves(s).bio) & isfinite(curves(s).tyr)); continue; end
    c = strainColor(cmap, curves(s).name);
    lw = 2.8; if strcmpi(curves(s).name, 'By4743 wt'); lw = 3.4; end
    plot(ax, curves(s).bio, curves(s).tyr, '-', 'LineWidth', lw, 'Color', c, 'DisplayName', curves(s).name);
    leg{end+1} = curves(s).name; %#ok<AGROW>
    maxB = max(maxB, max(curves(s).bio(isfinite(curves(s).bio))));
    maxT = max(maxT, max(curves(s).tyr(isfinite(curves(s).tyr))));
    if ~isempty(exp)
        j = find(strcmp({exp.name}, curves(s).name), 1);
        if ~isempty(j)
            scatter(ax, exp(j).bio, exp(j).tyr, 72, 'o', ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            maxB = max(maxB, max(exp(j).bio)); maxT = max(maxT, max(exp(j).tyr));
        end
    end
end
legend(ax, leg, 'Location', 'northeast', 'Interpreter', 'none', 'Box', 'off', 'TextColor', 'k', 'FontSize', 11);
xlabel(ax, 'Biomass yield [gDW/g glucose]', 'Color', 'k');
ylabel(ax, 'Tyrosol yield [g/g glucose]', 'Color', 'k');
xlim(ax, [0, max(0.01, 1.05 * maxB)]); ylim(ax, [0, max(0.001, 1.10 * maxT)]);
hold(ax, 'off'); box(ax, 'on');
savefig(fig, [out '.fig']);
exportgraphics(fig, [out '.png'], 'Resolution', 600, 'BackgroundColor', 'white');
exportgraphics(fig, [out '.svg'], 'ContentType', 'vector', 'BackgroundColor', 'white');
close(fig);
fprintf('  saved %s.{fig,png,svg}\n', out);
end

function exp = readFlaskYields(file,G)
T=readtable(file,'VariableNamingRule','preserve');
sc=col(T,{'STRAIN','Strain','strain'});
exp=repmat(struct('name','','bio',[],'tyr',[]),height(T),1);
for i=1:height(T)
    o=[T.(col(T,{'OD1'}))(i) T.(col(T,{'OD2'}))(i) T.(col(T,{'OD3'}))(i)];
    t=[T.(col(T,{'Tyrosol-1'}))(i) T.(col(T,{'Tyrosol-2'}))(i) T.(col(T,{'Tyrosol-3'}))(i)];
    exp(i).name=strtrim(char(string(T.(sc)(i))));
    exp(i).bio=(o+0.099)/3.6877/G; exp(i).tyr=t/(1000*G);
end
end

function writeExpCsv(exp, f)
strain = strings(0); rep = strings(0); bio = []; tyr = [];
for i = 1:numel(exp)
    for k = 1:numel(exp(i).bio)
        strain(end+1, 1) = string(exp(i).name); %#ok<AGROW>
        rep(end+1, 1) = "R" + k; %#ok<AGROW>
        bio(end+1, 1) = exp(i).bio(k); %#ok<AGROW>
        tyr(end+1, 1) = exp(i).tyr(k); %#ok<AGROW>
    end
end
writetable(table(strain, rep, bio, tyr, 'VariableNames', {'strain','rep','bioYield','tyrYield'}), f);
fprintf('  saved %s\n', f);
end

function writeCurveCsv(curves, f)
strain = strings(0); pt = []; bio = []; tyr = [];
for i = 1:numel(curves)
    n = numel(curves(i).bio);
    strain = [strain; repmat(string(curves(i).name), n, 1)]; %#ok<AGROW>
    pt = [pt; (1:n)']; %#ok<AGROW>
    bio = [bio; curves(i).bio(:)]; %#ok<AGROW>
    tyr = [tyr; curves(i).tyr(:)]; %#ok<AGROW>
end
writetable(table(strain, pt, bio, tyr, 'VariableNames', {'strain','point','bioYield','tyrYield'}), f);
fprintf('  saved %s\n', f);
end

function c = strainColors()
% gem12 (12-color tableau); wt = black, ByTOHdef = gray.
% ByTOH1-10 use gem12 slots 1-7 and 9-11 (skip slot 8 — it is gray).
gem12 = [0.1216 0.4667 0.7059; 0.6824 0.7804 0.9098; 0.5961 0.8745 0.5412; ...
    0.8392 0.1529 0.1569; 0.5804 0.4039 0.7412; 0.5490 0.3373 0.2941; ...
    0.8902 0.4667 0.7608; 0.4980 0.4980 0.4980; 0.7373 0.7412 0.1333; ...
    0.0902 0.7451 0.8118; 0.4157 0.2392 0.6039; 0.6941 0.3490 0.1569];
toh = {'ByTOH1','ByTOH2','ByTOH3','ByTOH4','ByTOH5','ByTOH6','ByTOH7','ByTOH8','ByTOH9','ByTOH10'};
tohGemIdx = [1 2 3 4 5 6 7 9 10 11];
c = containers.Map('KeyType', 'char', 'ValueType', 'any');
c('By4743 wt') = [0 0 0];
c('ByTOHdef') = [0.45 0.45 0.45];
for i = 1:numel(toh)
    c(toh{i}) = gem12(tohGemIdx(i), :);
end
end

function rgb = strainColor(cmap, strainName)
if isempty(cmap)
    rgb = [0.12 0.20 0.65];
    return;
end
name = strtrim(char(strainName));
if isKey(cmap, name)
    rgb = cmap(name);
else
    warning('No color for strain "%s" — using black.', name);
    rgb = [0 0 0];
end
end

function curves = loadCurveCsv(f)
T = readtable(f);
curves = repmat(struct('name','','bio',[],'tyr',[]), 0, 1);
for name = unique(T.strain, 'stable')'
    rows = strcmp(string(T.strain), string(name));
    curves(end+1).name = char(name); %#ok<AGROW>
    curves(end).bio = T.bioYield(rows);
    curves(end).tyr = T.tyrYield(rows);
end
end

function f=flaskFile()
h=char(java.lang.System.getProperty('user.home'));
p={fullfile(h,'Library','CloudStorage','OneDrive-Chalmers','Documents','tyrosol_ecYeasy','6. Cepas tirosol - rendimientos en matraz.xlsx'), ...
   fullfile(h,'Documents','tyrosol_ecYeasy','6. Cepas tirosol - rendimientos en matraz.xlsx')};
for i=1:2; if exist(p{i},'file'); f=p{i}; return; end; end; error('Flask workbook not found.');
end

function c=col(T,names)
for i=1:numel(names); h=T.Properties.VariableNames(strcmpi(T.Properties.VariableNames,names{i})); if ~isempty(h); c=h{1}; return; end; end
error('Column not found.');
end
