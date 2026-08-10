function mut = applyTyrosolMods(base, mods, umap, gur, oeCeil, kdScale)
% Apply OE/KD/KO edits. Used by envelopes and dFBA.
% oeCeil=true -> OE only raises ub (no forced floor). kdScale loosens KD caps.
if nargin < 4 || isempty(gur); gur = 1; end
if nargin < 5 || isempty(oeCeil); oeCeil = false; end
if nargin < 6 || isempty(kdScale); kdScale = 1; end

mut = base;
gi   = find(strcmpi(mut.rxnNames,'growth'),1);
gluc = find(strcmpi(mut.rxnNames,'D-glucose exchange (reversible)'),1);

for k = 1:size(mods,1)
    gene = char(mods{k,1});
    act  = char(mods{k,2});
    fac  = mods{k,3}; if iscell(fac); fac = fac{1}; end

    % already wired as pathway arms
    if strcmpi(act,'OE') && strcmpi(gene,'YDR380W') && any(strcmp(mut.rxns,'new_aro10_HPP')); continue; end
    if strcmpi(act,'OE') && strcmpi(gene,'YCR105W') && any(strcmp(mut.rxns,'new_adh7_tyrosol')); continue; end

    bu = 1e-9;
    if strcmpi(act,'OE') || strcmpi(act,'KD')
        if isKey(umap,gene); bu = umap(gene); else; bu = getUsage(mut,gene); end
    end

    prev = mut;
    mut = setGene(mut, gene, act, double(fac), bu, oeCeil, kdScale);

    if oeCeil && strcmpi(act,'OE'); continue; end

    if strcmpi(act,'KD') && ~canGrow(mut, gi, gluc, gur)
        nat = getUsageAt(prev, gene, gi, gluc, gur);
        if nat > 0
            mut = setGene(prev, gene, 'KD', 1.0, nat, oeCeil, kdScale);
            fprintf('      relax KD %s -> %.3g\n', gene, nat);
        end
        if ~canGrow(mut, gi, gluc, gur); mut = prev; end
    elseif strcmpi(act,'OE') && ~canGrow(mut, gi, gluc, gur)
        [mut, fr] = softenOE(prev, gene, bu, double(fac), gi, gluc, gur);
        fprintf('      relax OE %s -> %.2gx\n', gene, fr);
    end
end
end

function [mut, fr] = softenOE(prev, gene, bu, fac, gi, gluc, gur)
mut = prev; fr = 0;
ix = find(strcmpi(prev.enzGenes,gene),1);
if isempty(ix); return; end
rxs = find(contains(prev.rxnNames, prev.enzymes{ix}))';
for fr = [0.5 0.1 0]
    m = prev;
    for r = rxs
        if contains(m.rxnNames{r},'exchange') || m.ub(r) < 100
            m.ub(r) = bu*fac; m.lb(r) = 0;
        else
            m.lb(r) = fr*min(1000,bu*fac); m.ub(r) = 1000;
            if m.ub(r) <= m.lb(r); m.lb(r) = 0.99*m.ub(r); end
        end
    end
    if canGrow(m, gi, gluc, gur); mut = m; return; end
end
mut = m; fr = 0;
end

function ok = canGrow(m, gi, gluc, gur)
m = setParam(m,'ub',gluc,1.000001*gur);
m = setParam(m,'lb',gluc,0.999999*gur);
m = setParam(m,'obj',m.rxns{gi},1);
s = solveLP(m);
ok = ~isempty(s) && isfield(s,'x') && ~isempty(s.x) && s.x(gi) > 1e-9;
end

function u = getUsageAt(m, gene, gi, gluc, gur)
u = 0;
ix = find(strcmpi(m.enzGenes,gene),1); if isempty(ix); return; end
r = find(contains(m.rxnNames,m.enzymes{ix}),1); if isempty(r); return; end
m = setParam(m,'ub',gluc,1.000001*gur);
m = setParam(m,'lb',gluc,0.999999*gur);
s = solveLP(setParam(m,'obj',m.rxns{gi},1));
if ~isempty(s) && isfield(s,'x') && numel(s.x) >= r; u = max(s.x(r),0); end
end

function bu = getUsage(m, gene)
bu = 1e-9;
ix = find(strcmpi(m.enzGenes,gene),1); if isempty(ix); return; end
rx = find(contains(m.rxnNames,m.enzymes{ix}),1);
gi = find(strcmpi(m.rxnNames,'growth'),1);
if isempty(rx) || isempty(gi) || rx > numel(m.rxns); return; end
s = solveLP(setParam(m,'obj',m.rxns{gi},1));
if ~isempty(s) && isfield(s,'x') && numel(s.x) >= rx; bu = max(s.x(rx),1e-9); end
end

function m = setGene(m, gene, act, fac, bu, oeCeil, kdScale)
if nargin < 6 || isempty(oeCeil); oeCeil = false; end
if nargin < 7 || isempty(kdScale); kdScale = 1; end
ix = find(strcmpi(m.enzGenes,gene),1);
kf = 1; if strcmpi(act,'KD'); kf = kdScale; end

if strcmpi(act,'KO'); m = removeGenes(m,gene,false,false,false); end

if isempty(ix) && strcmpi(act,'KD')
    rx = find(contains(m.grRules,gene));
    if ~isempty(rx); m.S(:,rx) = m.S(:,rx)*fac; end
    return;
end
if isempty(ix); return; end

for r = find(contains(m.rxnNames, m.enzymes{ix}))'
    if contains(m.rxnNames{r},'exchange') || m.ub(r) < 100
        m.ub(r) = bu*fac*kf; m.lb(r) = 0;
    elseif fac > 1 && oeCeil
        m.ub(r) = 1000; m.lb(r) = 0;
    elseif fac > 1
        m.lb(r) = min(1000,bu*fac); m.ub(r) = 1000;
    else
        m.ub(r) = bu*fac*kf; m.lb(r) = 0;
    end
    if m.ub(r) <= m.lb(r); m.lb(r) = 0.99*m.ub(r); end
end
end
