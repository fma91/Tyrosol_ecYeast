%% Tyrosol predictive overflow dFBA — rho0 / no-mito ecModel batch cultures
%
%  addpath('scripts'); tyrosol_dfba
%
%  Forward (predictive) dynamic FBA for a batch culture on the rho0 / no-mito
%  ecModel. Only the initial state and two calibration targets come from the
%  experiment; glucose, biomass, ethanol and tyrosol trajectories are predicted.
%
%  Two batches are selectable through `experiment`:
%    * 'G5'       — ByTOHdef engineered strain (ARO1 OE, cumulative series)
%    * 'parental' — By4743 wt (petite) parental strain
%
%  Genotypes come from defineByTyrOHStrainTable + applyTyrosolMods (same
%  cumulative OE/KD/KO logic as the production envelopes, with OE floors
%  enforced — consistent with the envelope simulations).
%
%  Two parameters are auto-calibrated against the experimental batch:
%    1. KD looseness (kdScale) — matched to the measured post-lag growth rate.
%    2. Glucose uptake capacity (qSmax) — scanned to match the measured final
%       glucose (non-crash branch).
%
%  Spurious sugar/phosphate boundary exports are closed so glucose carbon
%  cannot leave the cell unmetabolised (closeGlucoseExport).

medium     = 'Min';
experiment = 'G5';   % 'G5' (ByTOHdef) | 'parental' (By4743 wt)

switch lower(experiment)
    case 'g5'
        strainName = 'ByTOHdef';
        expCsv     = 'experimental_G5_batch.csv';
        tStop      = 93;
        outTag     = 'G5';
        pushTyr    = true;   % engineered strain: maximise tyrosol within growth slack
        statCap    = inf;    % no explicit stationary-phase biomass cap
    case 'parental'
        strainName = 'By4743 wt';
        expCsv     = 'experimental_parental_batch.csv';
        tStop      = 72;
        outTag     = 'parental';
        pushTyr    = false;  % unengineered WT: no tyrosol pull
        statCap    = nan;    % cap growth at observed biomass plateau
    otherwise
        error('Unknown experiment: %s (use ''G5'' or ''parental'')', experiment);
end

dtGrid    = 1;              % h — coarse grid for seeding the uptake estimate
glcMW     = 0.180156;       % g/mmol
tyrMW     = 138.164/1000;   % g/mmol
ethMW     = 46.068/1000;    % g/mmol
vGlcUbMax = 1000;           % ceiling for the KD-relax sizing probe
KmGlc     = 0.5;            % g/L — Monod glucose half-saturation
growthDev = 0.10;           % growth-rate slack within which tyrosol is maximised
lagFrac   = 0.25;           % fraction of uptake capacity active during the lag phase
dtFwd     = 0.25;           % h — Euler integration step

pkg = fileparts(fileparts(mfilename('fullpath')));
addpath(pkg); addpath(fullfile(pkg,'scripts'));
home = char(java.lang.System.getProperty('user.home'));
addpath(genpath(fullfile(home,'Documents','CellFactory-ecYeastGEM','code')));
addpath(genpath(fullfile(home,'Documents','ecFactory','code')));
addpath(genpath(fullfile(home,'RAVEN')));
try; setRavenSolver('gurobi'); clear global RAVENSOLVER; catch; warning('gurobi unavailable; falling back to glpk'); end
set(0,'DefaultFigureVisible','off');

m   = loadModel(fullfile(pkg,'model','ecTyrosol.mat'));
m   = noMito(m);
tid = find(strcmpi(m.rxns,'new_tyrosol_ex'),1);
gi  = find(strcmpi(m.rxnNames,'growth'),1);
glc = find(strcmpi(m.rxnNames,'D-glucose exchange (reversible)'),1);
eth = find(strcmpi(m.rxns,'r_1761'),1);
umap = usageMap(fullfile(pkg,'results','candidates_L2.txt'));

dataDir = fullfile(pkg,'results','data');
figDir  = fullfile(pkg,'results','figures','dfba_no_mito');
if ~exist(dataDir,'dir'); mkdir(dataDir); end
if ~exist(figDir,'dir');  mkdir(figDir);  end

base    = baseModel(m, medium);
strains = defineByTyrOHStrainTable('measured', flaskFile());
strains = strains(ismember({strains.name}, {strainName}));
if isempty(strains); error('Strain %s not found.', strainName); end

% Size the KD-relax probe at the enzyme-limited uptake of the unconstrained model.
mProbe = setParam(base,'ub',glc,vGlcUbMax); mProbe = setParam(mProbe,'lb',glc,0);
sRef   = solveLP(setParam(mProbe,'obj',mProbe.rxns{gi},1), 0);
vRef   = 1;
if ~isempty(sRef) && isfield(sRef,'x') && ~isempty(sRef.x)
    vRef = max(abs(sRef.x(glc)), 1);
end

% Load and interpolate experimental batch data.
exp   = readExpBatch(fullfile(dataDir, expCsv));
tEnd  = min([tStop max(exp.glc.t) max(exp.bio.t)]);
tgrid = (0:dtGrid:tEnd)';
Gi    = interp1(exp.glc.t, exp.glc.v, tgrid, 'pchip');
Xi    = interp1(exp.bio.t, exp.bio.v, tgrid, 'pchip');

% Stationary-phase biomass cap (WT petite arrests; inf disables for G5).
if isnan(statCap)
    statCap = max(exp.bio.v(exp.bio.t <= tEnd + 1e-9));
end

% Product offsets: ethanol and tyrosol at t=0 from pre-culture carry-over.
eth0 = interp1(exp.eth.t, exp.eth.v, 0, 'pchip');
tyr0 = interp1(exp.tyr.t, exp.tyr.v, 0, 'pchip');

G0   = interp1(exp.glc.t, exp.glc.v, 0, 'pchip');
X0   = interp1(exp.bio.t, exp.bio.v, 0, 'pchip');
Xend = interp1(exp.bio.t, exp.bio.v, tEnd, 'pchip');
tgF  = (0:dtFwd:tEnd)';

% Lag: first time biomass departs the inoculum plateau.
iLag = find(exp.bio.v > 1.5*X0 & exp.bio.t <= tEnd + 1e-9, 1);
tLag = 0; if ~isempty(iLag); tLag = exp.bio.t(iLag); end

% Baseline specific uptake from the carbon balance (whole-batch average).
intX  = trapz(tgrid, Xi);
Gcons = max(Gi(1) - Gi(end), 1e-9);
qBase = (Gcons/glcMW) / intX;   % mmol glc/gDCW/h

%% --- Stage 1: calibrate KD looseness (kdScale) to match post-lag growth ---
muNeed = log(max(Xend, X0+1e-9)/X0) / max(tEnd-tLag, 1);
kdGrid = [1 2 3 4 5 6 7 8 10 12];
muCap  = zeros(numel(kdGrid),1);
fprintf('\n=== dFBA | %s | no_mito | ARO1 OE ===\n', strainName);
fprintf('  calibrating kdScale (target mu=%.4f 1/h)\n', muNeed);
for j = 1:numel(kdGrid)
    mj = applyTyrosolMods(base, strains(1).mods, umap, vRef, true, kdGrid(j));
    mp = setParam(mj,'lb',glc,0.999*qBase); mp = setParam(mp,'ub',glc,qBase);
    sj = solveLP(setParam(mp,'obj',mp.rxns{gi},1), 0);
    if okGrow(sj,gi); muCap(j) = sj.x(gi); end
    fprintf('  kdScale=%-3g -> mu=%.4f\n', kdGrid(j), muCap(j));
end
[~,ok] = min(abs(muCap - muNeed));
kdScale = kdGrid(ok); muRun = muCap(ok);
mutP = applyTyrosolMods(base, strains(1).mods, umap, vRef, true, kdScale);
fprintf('  -> kdScale=%g (mu=%.4f vs target %.4f)  lag=%.1f h\n', kdScale, muRun, muNeed, tLag);

%% --- Stage 2: calibrate qSmax to match measured final glucose ---
Gtarget  = interp1(exp.glc.t, exp.glc.v, tEnd, 'pchip');
tgCal    = (0:1.0:tEnd)';
boostGrid = 1.0:0.1:2.0;
fprintf('  scanning boost for glucose(%g h) closest to %.2f g/L  (qBase=%.2f)\n', tEnd, Gtarget, qBase);
bestErr = inf; uptakeBoost = boostGrid(1);
for b = boostGrid
    tr  = runOverflow(mutP, gi, glc, eth, tid, tgCal, G0, X0, glcMW, ethMW, tyrMW, ...
                      b*qBase, KmGlc, growthDev, tLag, lagFrac, pushTyr, statCap);
    gEnd = tr.glc(end);
    crash = gEnd < 1;
    fprintf('    boost=%.2f (qSmax=%.2f) -> glc(%g h)=%.2f%s\n', b, b*qBase, tEnd, gEnd, crash*' *CRASH*');
    if ~crash && abs(gEnd-Gtarget) < bestErr
        bestErr = abs(gEnd-Gtarget); uptakeBoost = b;
    end
end
qSmax = uptakeBoost * qBase;
fprintf('  -> boost=%.2f  qSmax=%.2f mmol/gDCW/h\n', uptakeBoost, qSmax);

%% --- Final run on fine grid ---
pred      = runOverflow(mutP, gi, glc, eth, tid, tgF, G0, X0, glcMW, ethMW, tyrMW, ...
                        qSmax, KmGlc, growthDev, tLag, lagFrac, pushTyr, statCap);
pred.name = strainName;
pred.eth  = pred.eth + eth0;
pred.tyr  = pred.tyr + tyr0;

fprintf('  [final] at %g h:  glc=%.2f (exp %.2f)  bio=%.3f (exp %.3f)  eth=%.2f (exp %.2f)  tyr=%.3f (exp %.3f)\n', ...
    tEnd, pred.glc(end), Gtarget, pred.bio(end), Xend, ...
    pred.eth(end), interp1(exp.eth.t, exp.eth.v, tEnd, 'linear','extrap'), ...
    pred.tyr(end), interp1(exp.tyr.t, exp.tyr.v, tEnd, 'linear','extrap'));

outBase = fullfile(figDir, sprintf('dfba_%s', outTag));
writeRunCsv(pred, pred.mu, pred.qS, fullfile(dataDir, sprintf('dfba_%s.csv', outTag)));
plotDFBA(pred, exp, outBase, 'Model (predictive dFBA)', ...
    sprintf('Predictive overflow dFBA — %s — no-mito', strainName));
fprintf('\nDone.\n');

%% =========================================================================
function ok = okGrow(s, gi)
ok = ~isempty(s) && isfield(s,'x') && ~isempty(s.x) && s.x(gi) > 1e-8;
end

%% --- overflow dFBA core ---
function traj = runOverflow(mut, gi, glc, eth, tid, tgrid, G0, X0, ...
                             glcMW, ethMW, tyrMW, qSmax, Km, growthDev, ...
                             tLag, lagFrac, pushTyr, statCap)
% Glucose uptake is forced to Monod kinetics; the enzyme budget spills excess
% carbon to ethanol (Crabtree overflow). Biomass, ethanol and tyrosol are
% predicted. During the lag phase the cell ferments at lagFrac of capacity
% with no net growth. In stationary phase growth arrests but fermentation
% continues. pushTyr maximises tyrosol within a 10% growth-rate slack.
if nargin < 15 || isempty(tLag);    tLag    = 0;    end
if nargin < 16 || isempty(lagFrac); lagFrac = 0;    end
if nargin < 17 || isempty(pushTyr); pushTyr = true; end
if nargin < 18 || isempty(statCap); statCap = inf;  end
n = numel(tgrid);
G = zeros(n,1); X = zeros(n,1); Eth = zeros(n,1); Tyr = zeros(n,1);
mu = zeros(n,1); qSu = zeros(n,1);
G(1) = G0; X(1) = X0; Cet = 0; Cty = 0;
t0 = tic;
for k = 1:n
    muK = 0; qSk = 0; vEt = 0; vTy = 0;
    Gk     = max(G(k), 0);
    inLag  = tgrid(k) < tLag;
    inStat = X(k) >= statCap;
    if inLag && lagFrac > 0
        qT = lagFrac * qSmax * Gk/(Km+Gk);
        mm = setParam(mut,'ub',glc,qT); mm = setParam(mm,'lb',glc,max(0,0.999*qT));
        mm = setParam(mm,'ub',mm.rxns{gi},1e-4);
        objLag = tid; if ~pushTyr; objLag = eth; end
        st = solveLP(setParam(mm,'obj',mm.rxns{objLag},1), 0);
        if ~isempty(st) && isfield(st,'x') && ~isempty(st.x) && st.stat > 0
            qSk = abs(st.x(glc)); vEt = max(0,st.x(eth)); vTy = max(0,st.x(tid));
        end
    elseif inStat
        qT = qSmax * Gk/(Km+Gk);
        mm = setParam(mut,'ub',glc,qT); mm = setParam(mm,'lb',glc,max(0,0.999*qT));
        mm = setParam(mm,'ub',mm.rxns{gi},1e-4);
        st = solveLP(setParam(mm,'obj',mm.rxns{eth},1), 0);
        if ~isempty(st) && isfield(st,'x') && ~isempty(st.x) && st.stat > 0
            qSk = abs(st.x(glc)); vEt = max(0,st.x(eth)); vTy = max(0,st.x(tid));
        end
    else
        qT = qSmax * Gk/(Km+Gk);
        mm = setParam(mut,'ub',glc,qT); mm = setParam(mm,'lb',glc,max(0,0.999*qT));
        sg = solveLP(setParam(mm,'obj',mm.rxns{gi},1), 0);
        if ~okGrow(sg,gi)   % forced uptake infeasible: fall back to cap-only
            mm = setParam(mut,'lb',glc,0); mm = setParam(mm,'ub',glc,qT);
            sg = solveLP(setParam(mm,'obj',mm.rxns{gi},1), 0);
        end
        if okGrow(sg,gi)
            muK = sg.x(gi); qSk = abs(sg.x(glc));
            if pushTyr
                vEt = max(0, sg.x(eth));
                mm2 = setParam(mm,'lb',mm.rxns{gi}, max(0,(1-growthDev)*muK));
                st  = solveLP(setParam(mm2,'obj',mm2.rxns{tid},1), 0);
                if ~isempty(st) && isfield(st,'x') && ~isempty(st.x) && st.stat > 0
                    vTy = max(0, st.x(tid));
                end
            else
                sp = solveLP(setParam(mm,'obj',mm.rxns{gi},1), 1);
                if okGrow(sp,gi); vEt = max(0,sp.x(eth)); vTy = max(0,sp.x(tid));
                else;             vEt = max(0,sg.x(eth)); end
            end
        end
    end
    mu(k) = muK; qSu(k) = qSk;
    if k < n
        dt     = tgrid(k+1) - tgrid(k);
        X(k+1) = X(k) + muK*X(k)*dt;
        G(k+1) = max(0, G(k) - qSk*X(k)*glcMW*dt);
        Cet    = Cet + vEt*X(k)*dt;
        Cty    = Cty + vTy*X(k)*dt;
    end
    Eth(k) = Cet*ethMW; Tyr(k) = Cty*tyrMW;
    if mod(k-1,20)==0 || k==n
        fprintf('    t=%5.1f h  G=%.2f  X=%.3f  mu=%.4f  qS=%.2f  vEtOH=%.3g  (%.1fs)\n', ...
            tgrid(k), G(k), X(k), muK, qSk, vEt, toc(t0));
    end
end
traj.t = tgrid; traj.glc = G; traj.bio = X;
traj.eth = Eth; traj.tyr = Tyr; traj.mu = mu; traj.qS = qSu;
end

%% --- plotting ---
function plotDFBA(traj, exp, out, modelLabel, titleStr)
if nargin < 4 || isempty(modelLabel); modelLabel = 'Model'; end
if nargin < 5; titleStr = ''; end
vars  = {'glc','bio','eth','tyr'};
expv  = {exp.glc, exp.bio, exp.eth, exp.tyr};
ylabs = {'Glucose [g/L]','Biomass [gDCW/L]','Ethanol [g/L]','Tyrosol [g/L]'};
have  = isfield(traj, vars);
vars = vars(have); expv = expv(have); ylabs = ylabs(have);
np = numel(vars);
modelCol = [0.12 0.20 0.65];
tMax = max(traj.t);
fig = figure('Color','w','Visible','off','Position',[100 100 480*np 480]);
for p = 1:np
    ax = subplot(1,np,p); hold(ax,'on');
    set(ax,'Color','w','XColor','k','YColor','k','LineWidth',1.4,'FontSize',12, ...
        'TickDir','out','TickLength',[0.018 0.025]);
    hM = plot(ax, traj.t, traj.(vars{p}), '-', 'LineWidth',2.8, 'Color',modelCol);
    e    = expv{p};
    keep = e.t <= tMax + 1e-9;
    hE = errorbar(ax, e.t(keep), e.v(keep), e.sd(keep), 'o', 'MarkerSize',6, ...
        'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','k', ...
        'Color',[0.85 0.33 0.10], 'LineWidth',1.1, 'CapSize',4);
    xlabel(ax,'Time [h]','Color','k'); ylabel(ax,ylabs{p},'Color','k');
    xlim(ax,[0 tMax]); box(ax,'on'); hold(ax,'off');
    if p == 1
        legend(ax, [hM hE], {modelLabel,'Experiment'}, ...
            'Location','northeast','Box','off','TextColor','k','FontSize',11);
    end
end
if ~isempty(titleStr)
    sgtitle(sprintf('%s  —  %s', traj.name, titleStr), 'FontSize',14,'FontWeight','bold');
end
savefig(fig,[out '.fig']);
exportgraphics(fig,[out '.png'],'Resolution',600,'BackgroundColor','white');
exportgraphics(fig,[out '.svg'],'ContentType','vector','BackgroundColor','white');
close(fig);
fprintf('  saved %s.{fig,png,svg}\n', out);
end

%% --- I/O helpers ---
function exp = readExpBatch(f)
T = readtable(f);
exp.glc = pickVar(T,'glucose');
exp.bio = pickVar(T,'biomass');
exp.tyr = pickVar(T,'tyrosol');
exp.eth = pickVar(T,'ethanol');
end

function s = pickVar(T, name)
rows = strcmp(string(T.variable), name);
s.t  = T.time_h(rows); s.v = T.value_gL(rows); s.sd = T.sd_gL(rows);
[s.t,ord] = sort(s.t); s.v = s.v(ord); s.sd = s.sd(ord);
end

function writeRunCsv(traj, muV, qsV, f)
eth = zeros(size(traj.t)); if isfield(traj,'eth'); eth = traj.eth; end
writetable(table(traj.t, traj.glc, traj.bio, eth, traj.tyr, muV(:), qsV(:), ...
    'VariableNames',{'time_h','glucose_gL','biomass_gL','ethanol_gL','tyrosol_gL','mu_1h','qS_mmol_gDCW_h'}), f);
fprintf('  saved %s\n', f);
end

%% --- model helpers ---
function m = loadModel(f)
raw = load(f); fn = fieldnames(raw); m = raw.(fn{1});
for fld = {'lb','ub','c','b','rev'}
    if isfield(m,fld{1}); m.(fld{1}) = full(double(m.(fld{1})(:))); end
end
for fld = {'rxns','rxnNames','genes','grRules'}
    if isfield(m,fld{1})
        v = m.(fld{1}); if isstring(v)||ischar(v); v=cellstr(v); end
        m.(fld{1}) = v(:);
    end
end
if isfield(m,'S') && ~issparse(m.S); m.S = sparse(m.S); end
end

function b = baseModel(m, medium)
b = changeMedia_batch(m,'D-glucose exchange (reversible)',medium);
b = setParam(b,'lb',find(strcmpi(b.rxnNames,'growth')),0);
if any(strcmpi(b.rxns,'r_2111')); b=setParam(b,'lb','r_2111',0); b=setParam(b,'ub','r_2111',1000); end
b = closeGlucoseExport(b);
end

function b = closeGlucoseExport(b)
% Block secretion through spurious sugar/phosphate boundary exchanges so that
% glucose carbon cannot leave the cell unmetabolised. Glucose uptake and
% inorganic phosphate are left untouched.
leak = {'r_4502','r_4504','r_4507','r_4538','r_4543','r_4539','r_4547', ...
        'r_1651','r_4499','r_4522','r_4535','r_1709','r_1715','r_1716','r_1650'};
for i = 1:numel(leak)
    if any(strcmp(b.rxns,leak{i})); b = setParam(b,'ub',leak{i},0); end
end
end

function m = noMito(m)
for g = {'Q0045','Q0080','Q0085','Q0105','Q0130','Q0250','Q0275'}
    if any(strcmpi(m.genes,g{1}))||any(strcmpi(m.enzGenes,g{1}))
        m = removeGenes(m,g{1},false,false,false);
    end
end
end

function u = usageMap(f)
u = containers.Map; if ~exist(f,'file'); return; end
T = readtable(f,'FileType','text','Delimiter','\t');
for i=1:height(T)
    v=T.maxUsageBio(i); if isnan(v)||v<=0; v=max(T.maxUsage(i),1e-9); end
    u(char(T.genes(i)))=1.01*v;
end
end

function f = flaskFile()
h = char(java.lang.System.getProperty('user.home'));
p = {fullfile(h,'Library','CloudStorage','OneDrive-Chalmers','Documents','tyrosol_ecYeasy','6. Cepas tirosol - rendimientos en matraz.xlsx'), ...
     fullfile(h,'Documents','tyrosol_ecYeasy','6. Cepas tirosol - rendimientos en matraz.xlsx')};
for i = 1:2; if exist(p{i},'file'); f=p{i}; return; end; end
error('Flask workbook not found.');
end
