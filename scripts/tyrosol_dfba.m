%% Predictive overflow dFBA (no-mito)
% experiment = 'G5' (ByTOHdef) or 'parental' (By4743 wt)
% addpath('scripts'); tyrosol_dfba

medium     = 'Min';
experiment = 'G5';

switch lower(experiment)
    case 'g5'
        strain  = 'ByTOHdef';
        expFile = 'experimental_G5_batch.csv';
        tStop   = 93;
        tag     = 'G5';
        pushTyr = true;
        Xmax    = inf;
    case 'parental'
        strain  = 'By4743 wt';
        expFile = 'experimental_parental_batch.csv';
        tStop   = 72;
        tag     = 'parental';
        pushTyr = false;
        Xmax    = nan;          % use measured plateau
    otherwise
        error('Use ''G5'' or ''parental''.');
end

% Optional: extendTo continues the same calibrated trajectory past tStop
% (e.g. tyrosol_dfba_full sets extendTo=167 → fit on 0–93 h, run to 167 h).
tFit = tStop;
tRun = tStop;
if exist('extendTo', 'var') && ~isempty(extendTo)
    tRun = extendTo;
    tag  = [tag '_full'];
end

dt      = 0.25;             % h
glcMW   = 0.180156;
tyrMW   = 138.164/1000;
ethMW   = 46.068/1000;
Km      = 0.5;              % g/L
muSlack = 0.10;
lagFrac = 0.25;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root,'scripts'));
home = char(java.lang.System.getProperty('user.home'));
addpath(genpath(fullfile(home,'Documents','CellFactory-ecYeastGEM','code')));
addpath(genpath(fullfile(home,'Documents','ecFactory','code')));
addpath(genpath(fullfile(home,'RAVEN')));
try; setRavenSolver('gurobi'); clear global RAVENSOLVER;
catch; warning('Gurobi not found, using default solver.'); end
set(0,'DefaultFigureVisible','off');

model = loadEcModel(fullfile(root,'model','ecTyrosol.mat'));
model = knockOutMito(model);
tid = find(strcmpi(model.rxns,'new_tyrosol_ex'),1);
gi  = find(strcmpi(model.rxnNames,'growth'),1);
glc = find(strcmpi(model.rxnNames,'D-glucose exchange (reversible)'),1);
eth = find(strcmpi(model.rxns,'r_1761'),1);
umap = loadUsageMap(fullfile(root,'results','candidates_L2.txt'));

dataDir = fullfile(root,'results','data');
figDir  = fullfile(root,'results','figures','dfba_no_mito');
if ~exist(dataDir,'dir'); mkdir(dataDir); end
if ~exist(figDir,'dir'); mkdir(figDir); end

base = setMinGlucose(model, medium);
base = closeLeaks(base);

strains = defineByTyrOHStrainTable('measured', findYieldsFile());
strains = strains(ismember({strains.name}, {strain}));
if isempty(strains); error('Strain %s not in table.', strain); end

% uptake capacity of unconstrained model (for KD sizing)
tmp = setParam(base,'ub',glc,1000); tmp = setParam(tmp,'lb',glc,0);
sol = solveLP(setParam(tmp,'obj',tmp.rxns{gi},1),0);
vRef = 1;
if ~isempty(sol) && isfield(sol,'x') && ~isempty(sol.x)
    vRef = max(abs(sol.x(glc)),1);
end

exp  = readBatch(fullfile(dataDir, expFile));
tFitEnd = min([tFit max(exp.glc.t) max(exp.bio.t)]);
tRunEnd = min([tRun max(exp.glc.t) max(exp.bio.t)]);
t1   = (0:1:tFitEnd)';
Gi   = interp1(exp.glc.t, exp.glc.v, t1, 'pchip');
Xi   = interp1(exp.bio.t, exp.bio.v, t1, 'pchip');

if isnan(Xmax); Xmax = max(exp.bio.v(exp.bio.t <= tFitEnd)); end

G0   = interp1(exp.glc.t, exp.glc.v, 0, 'pchip');
X0   = interp1(exp.bio.t, exp.bio.v, 0, 'pchip');
Xfit = interp1(exp.bio.t, exp.bio.v, tFitEnd, 'pchip');
Xend = interp1(exp.bio.t, exp.bio.v, tRunEnd, 'pchip');
eth0 = interp1(exp.eth.t, exp.eth.v, 0, 'pchip');
tyr0 = interp1(exp.tyr.t, exp.tyr.v, 0, 'pchip');
tF   = (0:dt:tRunEnd)';

% lag = first time biomass rises above inoculum (fit window)
iLag = find(exp.bio.v > 1.5*X0 & exp.bio.t <= tFitEnd, 1);
tLag = 0; if ~isempty(iLag); tLag = exp.bio.t(iLag); end

qBase = max(Gi(1)-Gi(end),1e-9)/glcMW / trapz(t1,Xi);   % mmol/gDCW/h

%% 1) fit kdScale to post-lag growth (fit horizon only)
muNeed = log(max(Xfit,X0+1e-9)/X0) / max(tFitEnd-tLag,1);
kdList = [1 2 3 4 5 6 7 8 10 12];
muList = zeros(size(kdList));
fprintf('\n=== dFBA %s (no-mito) ===\n', strain);
if tRunEnd > tFitEnd
    fprintf('  fit horizon: 0–%.0f h; integrate to %.0f h (same parameters)\n', tFitEnd, tRunEnd);
end
fprintf('  kdScale scan (target mu=%.4f)\n', muNeed);
for j = 1:numel(kdList)
    m = applyTyrosolMods(base, strains(1).mods, umap, vRef, true, kdList(j));
    m = setParam(m,'lb',glc,0.999*qBase); m = setParam(m,'ub',glc,qBase);
    s = solveLP(setParam(m,'obj',m.rxns{gi},1),0);
    if okSol(s,gi); muList(j) = s.x(gi); end
    fprintf('  kd=%-3g  mu=%.4f\n', kdList(j), muList(j));
end
[~,j] = min(abs(muList - muNeed));
kd = kdList(j);
mut = applyTyrosolMods(base, strains(1).mods, umap, vRef, true, kd);
fprintf('  -> kdScale=%g (mu=%.4f, need %.4f), lag=%.1f h\n', kd, muList(j), muNeed, tLag);

%% 2) fit uptake boost to glucose at fit horizon (NOT tRunEnd)
Gtarget = interp1(exp.glc.t, exp.glc.v, tFitEnd, 'pchip');
boosts  = 1.0:0.1:2.0;
bestErr = inf; bestB = boosts(1);
fprintf('  boost scan (target G(%.0f h)=%.2f, qBase=%.2f)\n', tFitEnd, Gtarget, qBase);
for b = boosts
    tr = runDFBA(mut, gi, glc, eth, tid, t1, G0, X0, glcMW, ethMW, tyrMW, ...
                 b*qBase, Km, muSlack, tLag, lagFrac, pushTyr, Xmax);
    crash = tr.glc(end) < 1;
    note = ''; if crash; note = ' *crash*'; end
    fprintf('    boost=%.2f  G=%.2f%s\n', b, tr.glc(end), note);
    if ~crash && abs(tr.glc(end)-Gtarget) < bestErr
        bestErr = abs(tr.glc(end)-Gtarget); bestB = b;
    end
end
qSmax = bestB * qBase;
fprintf('  -> boost=%.2f  qSmax=%.2f\n', bestB, qSmax);

%% final (may extend past fit horizon with locked parameters)
pred = runDFBA(mut, gi, glc, eth, tid, tF, G0, X0, glcMW, ethMW, tyrMW, ...
               qSmax, Km, muSlack, tLag, lagFrac, pushTyr, Xmax);
pred.name = strain;
pred.eth  = pred.eth + eth0;
pred.tyr  = pred.tyr + tyr0;

GexpEnd = interp1(exp.glc.t, exp.glc.v, tRunEnd, 'pchip');
fprintf('  @fit %.0fh: G=%.2f (exp %.2f)\n', tFitEnd, ...
    interp1(pred.t, pred.glc, tFitEnd, 'linear'), Gtarget);
fprintf('  final @%gh: G=%.2f (exp %.2f) X=%.3f (exp %.3f) EtOH=%.2f Tyr=%.3f\n', ...
    tRunEnd, pred.glc(end), GexpEnd, pred.bio(end), Xend, pred.eth(end), pred.tyr(end));

out = fullfile(figDir, sprintf('dfba_%s', tag));
writeDFBA(pred, fullfile(dataDir, sprintf('dfba_%s.csv', tag)));
plotDFBA(pred, exp, out);
fprintf('\nDone.\n');

%% ---
function ok = okSol(s, gi)
ok = ~isempty(s) && isfield(s,'x') && ~isempty(s.x) && s.x(gi) > 1e-8;
end

function tr = runDFBA(mut, gi, glc, eth, tid, t, G0, X0, glcMW, ethMW, tyrMW, ...
                      qSmax, Km, muSlack, tLag, lagFrac, pushTyr, Xmax)
if nargin < 15 || isempty(tLag); tLag = 0; end
if nargin < 16 || isempty(lagFrac); lagFrac = 0; end
if nargin < 17 || isempty(pushTyr); pushTyr = true; end
if nargin < 18 || isempty(Xmax); Xmax = inf; end

n = numel(t);
G = zeros(n,1); X = zeros(n,1); Et = zeros(n,1); Ty = zeros(n,1);
mu = zeros(n,1); qS = zeros(n,1);
G(1) = G0; X(1) = X0;
cEt = 0; cTy = 0;
t0 = tic;

for k = 1:n
    muK = 0; q = 0; vE = 0; vT = 0;
    Gk = max(G(k),0);
    lag  = t(k) < tLag;
    stop = X(k) >= Xmax;

    if lag && lagFrac > 0
        qT = lagFrac * qSmax * Gk/(Km+Gk);
        m = setParam(mut,'ub',glc,qT); m = setParam(m,'lb',glc,max(0,0.999*qT));
        m = setParam(m,'ub',m.rxns{gi},1e-4);
        obj = tid; if ~pushTyr; obj = eth; end
        s = solveLP(setParam(m,'obj',m.rxns{obj},1),0);
        if ~isempty(s) && isfield(s,'x') && ~isempty(s.x) && s.stat > 0
            q = abs(s.x(glc)); vE = max(0,s.x(eth)); vT = max(0,s.x(tid));
        end
    elseif stop
        qT = qSmax * Gk/(Km+Gk);
        m = setParam(mut,'ub',glc,qT); m = setParam(m,'lb',glc,max(0,0.999*qT));
        m = setParam(m,'ub',m.rxns{gi},1e-4);
        s = solveLP(setParam(m,'obj',m.rxns{eth},1),0);
        if ~isempty(s) && isfield(s,'x') && ~isempty(s.x) && s.stat > 0
            q = abs(s.x(glc)); vE = max(0,s.x(eth)); vT = max(0,s.x(tid));
        end
    else
        qT = qSmax * Gk/(Km+Gk);
        m = setParam(mut,'ub',glc,qT); m = setParam(m,'lb',glc,max(0,0.999*qT));
        sg = solveLP(setParam(m,'obj',m.rxns{gi},1),0);
        if ~okSol(sg,gi)
            m = setParam(mut,'lb',glc,0); m = setParam(m,'ub',glc,qT);
            sg = solveLP(setParam(m,'obj',m.rxns{gi},1),0);
        end
        if okSol(sg,gi)
            muK = sg.x(gi); q = abs(sg.x(glc));
            if pushTyr
                vE = max(0,sg.x(eth));
                m2 = setParam(m,'lb',m.rxns{gi}, max(0,(1-muSlack)*muK));
                st = solveLP(setParam(m2,'obj',m2.rxns{tid},1),0);
                if ~isempty(st) && isfield(st,'x') && ~isempty(st.x) && st.stat > 0
                    vT = max(0,st.x(tid));
                end
            else
                sp = solveLP(setParam(m,'obj',m.rxns{gi},1),1);
                if okSol(sp,gi); vE = max(0,sp.x(eth)); vT = max(0,sp.x(tid));
                else; vE = max(0,sg.x(eth)); end
            end
        end
    end

    mu(k) = muK; qS(k) = q;
    if k < n
        d = t(k+1)-t(k);
        X(k+1) = X(k) + muK*X(k)*d;
        G(k+1) = max(0, G(k) - q*X(k)*glcMW*d);
        cEt = cEt + vE*X(k)*d;
        cTy = cTy + vT*X(k)*d;
    end
    Et(k) = cEt*ethMW; Ty(k) = cTy*tyrMW;

    if mod(k-1,20)==0 || k==n
        fprintf('    t=%5.1f  G=%.2f  X=%.3f  mu=%.4f  (%.0fs)\n', t(k), G(k), X(k), muK, toc(t0));
    end
end

tr.t = t; tr.glc = G; tr.bio = X; tr.eth = Et; tr.tyr = Ty; tr.mu = mu; tr.qS = qS;
end

function plotDFBA(tr, exp, out)
vars = {'glc','bio','eth','tyr'};
labs = {'Glucose [g/L]','Biomass [gDCW/L]','Ethanol [g/L]','Tyrosol [g/L]'};
exps = {exp.glc, exp.bio, exp.eth, exp.tyr};
n = numel(vars);
fig = figure('Color','w','Visible','off','Position',[100 100 480*n 480]);
tMax = max(tr.t);
for p = 1:n
    ax = subplot(1,n,p); hold(ax,'on');
    set(ax,'LineWidth',1.4,'FontSize',12,'TickDir','out');
    hM = plot(ax, tr.t, tr.(vars{p}), '-', 'LineWidth',2.8, 'Color',[0.12 0.20 0.65]);
    e = exps{p}; keep = e.t <= tMax;
    hE = errorbar(ax, e.t(keep), e.v(keep), e.sd(keep), 'o', 'MarkerSize',6, ...
        'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','k', ...
        'Color',[0.85 0.33 0.10], 'LineWidth',1.1, 'CapSize',4);
    xlabel(ax,'Time [h]'); ylabel(ax,labs{p});
    xlim(ax,[0 tMax]); box(ax,'on'); hold(ax,'off');
    if p==1
        legend(ax,[hM hE],{'Model','Experiment'},'Location','northeast','Box','off');
    end
end
sgtitle(sprintf('%s — predictive dFBA (no-mito)', tr.name), 'FontSize',14, 'FontWeight','bold');
savefig(fig,[out '.fig']);
exportgraphics(fig,[out '.png'],'Resolution',600,'BackgroundColor','white');
exportgraphics(fig,[out '.svg'],'ContentType','vector','BackgroundColor','white');
close(fig);
fprintf('  saved %s.{fig,png,svg}\n', out);
end

function exp = readBatch(f)
T = readtable(f);
exp.glc = getVar(T,'glucose');
exp.bio = getVar(T,'biomass');
exp.tyr = getVar(T,'tyrosol');
exp.eth = getVar(T,'ethanol');
end

function s = getVar(T, name)
rows = strcmp(string(T.variable), name);
s.t = T.time_h(rows); s.v = T.value_gL(rows); s.sd = T.sd_gL(rows);
[s.t,ord] = sort(s.t); s.v = s.v(ord); s.sd = s.sd(ord);
end

function writeDFBA(tr, f)
eth = zeros(size(tr.t)); if isfield(tr,'eth'); eth = tr.eth; end
writetable(table(tr.t,tr.glc,tr.bio,eth,tr.tyr,tr.mu(:),tr.qS(:), ...
    'VariableNames',{'time_h','glucose_gL','biomass_gL','ethanol_gL','tyrosol_gL','mu_1h','qS_mmol_gDCW_h'}), f);
fprintf('  saved %s\n', f);
end

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
    b = setParam(b,'lb','r_2111',0); b = setParam(b,'ub','r_2111',1000);
end
end

function b = closeLeaks(b)
% block spurious sugar/phosphate exports
rx = {'r_4502','r_4504','r_4507','r_4538','r_4543','r_4539','r_4547', ...
      'r_1651','r_4499','r_4522','r_4535','r_1709','r_1715','r_1716','r_1650'};
for i = 1:numel(rx)
    if any(strcmp(b.rxns,rx{i})); b = setParam(b,'ub',rx{i},0); end
end
end

function m = knockOutMito(m)
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

function f = findYieldsFile()
h = char(java.lang.System.getProperty('user.home'));
c = {fullfile(h,'Library','CloudStorage','OneDrive-Chalmers','Documents','tyrosol_ecYeasy','6. Cepas tirosol - rendimientos en matraz.xlsx'), ...
     fullfile(h,'Documents','tyrosol_ecYeasy','6. Cepas tirosol - rendimientos en matraz.xlsx')};
for i = 1:2
    if exist(c{i},'file'); f = c{i}; return; end
end
error('Yields Excel not found.');
end
