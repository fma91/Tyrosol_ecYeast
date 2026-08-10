function mut = applyTyrosolMods(base, mods, umap, GUR, oeCeiling, kdScale)
% APPLYTYROSOLMODS  Apply a strain's cumulative OE/KD/KO edits to an ecModel.
%   Shared by tyrosol_envelopes.m and tyrosol_dfba.m so the genotype logic
%   stays identical between the steady-state envelopes and the dynamic runs.
%
%   OE  -> force a flux floor on the enzyme usage (base_usage*factor) with a
%          1000 ceiling; KD -> cap usage at base_usage*factor; KO -> removeGenes.
%   Any KD/OE that drives growth to ~0 on minimal medium is auto-relaxed:
%   KDs are reopened to the enzyme's natural usage, OE floors are lowered
%   (0.5x, 0.1x, ceiling-only) until the model grows again.
%
%   GUR (optional, default 1) is the glucose uptake the feasibility checks and
%   the essential-KD relaxation are sized at. Envelopes use 1 (their fixed
%   normalisation); the dFBA passes the model's enzyme-limited uptake so the
%   relaxed essential-KD caps don't strangle growth at high uptake.
if nargin < 4 || isempty(GUR); GUR = 1; end
% oeCeiling=true makes OE a ceiling-only edit (raise ub to 1000, no forced flux
% floor). Used by the data-constrained dFBA so OE edits don't impose minimum
% fluxes that make the LP degenerate/slow at the low measured glucose uptake.
% Default false preserves the envelope behaviour (forced OE floor + auto-relax).
if nargin < 5 || isempty(oeCeiling); oeCeiling = false; end
% kdScale (optional, default 1) multiplies every KD usage cap (bu*f*kdScale).
% >1 loosens the knockdowns so the engineered strain is not throttled below its
% measured growth in a free-running dFBA; 1 reproduces the original KD strength.
if nargin < 6 || isempty(kdScale); kdScale = 1; end
mut = base;
gi = find(strcmpi(mut.rxnNames,'growth'),1);
gluc = find(strcmpi(mut.rxnNames,'D-glucose exchange (reversible)'),1);
for k=1:size(mods,1)
    g=char(mods{k,1}); a=char(mods{k,2}); f=mods{k,3}; if iscell(f); f=f{1}; end
    if strcmpi(a,'OE')&&strcmpi(g,'YDR380W')&&any(strcmp(mut.rxns,'new_aro10_HPP')); continue; end
    if strcmpi(a,'OE')&&strcmpi(g,'YCR105W')&&any(strcmp(mut.rxns,'new_adh7_tyrosol')); continue; end
    bu=1e-9;
    if strcmpi(a,'OE') || strcmpi(a,'KD')
        if isKey(umap,g); bu=umap(g); else; bu=enzUsage(mut,g); end
    end
    prev = mut;
    mut = patchGene(mut,{g,a,double(f)},bu,oeCeiling,kdScale);
    if oeCeiling && strcmpi(a,'OE')
        continue;   % ceiling-only OE never forces a floor, so no relax needed
    end
    if strcmpi(a,'KD') && ~growsAt(mut, gi, gluc, GUR)
        nat = enzUsageAt(prev, g, gi, gluc, GUR);
        if nat > 0
            mut = patchGene(prev, {g,'KD',1.0}, nat, oeCeiling, kdScale);
            fprintf('      [auto-relax] KD %s -> cap at natural usage %.3g (kept feasible)\n', g, nat);
        end
        if ~growsAt(mut, gi, gluc, GUR)
            mut = prev;  % relaxation did not help: cause is elsewhere, drop this KD
        end
    elseif strcmpi(a,'OE') && ~growsAt(mut, gi, gluc, GUR)
        [mut, fr] = relaxOE(prev, g, bu, double(f), gi, gluc, GUR);
        fprintf('      [auto-relax] OE %s -> flux floor reduced to %.2gx (kept feasible)\n', g, fr);
    end
end
end

function [mut, frUsed] = relaxOE(prev, g, bu, f, gi, gluc, GUR)
% Reduce the forced OE lower-bound floor (down to ceiling-only) until the
% model grows again, keeping as much overexpression forcing as feasible.
mut = prev; frUsed = 0;
ix = find(strcmpi(prev.enzGenes,g),1);
if isempty(ix); return; end
rxs = find(contains(prev.rxnNames,prev.enzymes{ix}))';
for fr = [0.5 0.1 0]
    mm = prev;
    for r = rxs
        if contains(mm.rxnNames{r},'exchange') || mm.ub(r)<100
            mm.ub(r)=bu*f; mm.lb(r)=0;
        else
            mm.lb(r)=fr*min(1000,bu*f); mm.ub(r)=1000;
            if mm.ub(r)<=mm.lb(r); mm.lb(r)=0.99*mm.ub(r); end
        end
    end
    if growsAt(mm, gi, gluc, GUR); mut=mm; frUsed=fr; return; end
end
mut = mm; frUsed = 0;  % ceiling-only as last resort
end

function ok = growsAt(m, gi, gluc, GUR)
m = setParam(m,'ub',gluc,1.000001*GUR);
m = setParam(m,'lb',gluc,0.999999*GUR);
m = setParam(m,'obj',m.rxns{gi},1);
s = solveLP(m);
ok = ~isempty(s) && isfield(s,'x') && ~isempty(s.x) && s.x(gi) > 1e-9;
end

function u = enzUsageAt(m, g, gi, gluc, GUR)
u = 0;
ix = find(strcmpi(m.enzGenes,g),1); if isempty(ix); return; end
r = find(contains(m.rxnNames,m.enzymes{ix}),1); if isempty(r); return; end
m = setParam(m,'ub',gluc,1.000001*GUR);
m = setParam(m,'lb',gluc,0.999999*GUR);
m = setParam(m,'obj',m.rxns{gi},1);
s = solveLP(m);
if ~isempty(s) && isfield(s,'x') && numel(s.x)>=r; u = max(s.x(r), 0); end
end

function bu = enzUsage(m,g)
bu=1e-9; ix=find(strcmpi(m.enzGenes,g),1); if isempty(ix); return; end
rx=find(contains(m.rxnNames,m.enzymes{ix}),1);
gi=find(strcmpi(m.rxnNames,'growth'),1);
if isempty(rx) || isempty(gi) || rx > numel(m.rxns); return; end
sol=solveLP(setParam(m,'obj',m.rxns{gi},1));
if ~isempty(sol)&&isfield(sol,'x')&&numel(sol.x)>=rx; bu=max(sol.x(rx),1e-9); end
end

function m = patchGene(m,edit,bu,oeCeiling,kdScale)
if nargin < 4 || isempty(oeCeiling); oeCeiling = false; end
if nargin < 5 || isempty(kdScale); kdScale = 1; end
g=edit{1}; a=edit{2}; f=edit{3}; ix=find(strcmpi(m.enzGenes,g),1);
kf = 1; if strcmpi(a,'KD'); kf = kdScale; end   % only KD caps are loosened
if strcmpi(a,'KO'); m=removeGenes(m,g,false,false,false); end
if isempty(ix)&&strcmpi(a,'KD'); rx=find(contains(m.grRules,g)); if ~isempty(rx); m.S(:,rx)=m.S(:,rx)*f; end; return; end
if isempty(ix); return; end
for r=find(contains(m.rxnNames,m.enzymes{ix}))'
    if contains(m.rxnNames{r},'exchange')||m.ub(r)<100; m.ub(r)=bu*f*kf; m.lb(r)=0;
    elseif f>1 && oeCeiling; m.ub(r)=1000; m.lb(r)=0;            % ceiling-only OE
    elseif f>1; m.lb(r)=min(1000,bu*f); m.ub(r)=1000; else; m.ub(r)=bu*f*kf; m.lb(r)=0; end
    if m.ub(r)<=m.lb(r); m.lb(r)=0.99*m.ub(r); end
end
end
