"""Build ecTyrosol.mat — enzyme-constrained model for tyrosol production.

Starting point: ecYeastGEM_batch.mat (S. cerevisiae ecModel, CellFactory-ecYeastGEM).

Modeling assumptions (documented in code and in docs/METHODS.md):

  A1. Tyrosol is produced via the native Ehrlich pathway from tyrosine:
      L-tyrosine -> 4-hydroxyphenylpyruvate (4-HPP) -> 4-hydroxyphenylacetaldehyde
      (4-HPAA) -> tyrosol. This route uses endogenous aromatic metabolism rather
      than a heterologous tyrosine-to-tyramine bypass.

  A2. Decarboxylation of 4-HPP is assigned to ARO10 (YDR380W), consistent with
      reported promiscuity of Aro10p toward aromatic alpha-keto acids.

  A3. Reduction of 4-HPAA to tyrosol is assigned to ADH7 (YCR105W, NADPH-
      dependent aldehyde reductase), the native yeast enzyme class for this
      chemistry in related alcohol pathways.

  A4. Enzyme capacity on the new pathway steps is represented with GECKO-style
      arm reactions and kcat = 1000 1/s on ARO10 and ADH7 branches. This value
      is a modeling assumption for a dedicated, high-flux production route
      (S-coefficient = 1 / (kcat * 3600)).

  A5. Wild-type genetic background: no chassis gene deletions are applied in
      this model build.

  A6. Product export is modeled as cytoplasmic tyrosol, transport to the
      extracellular space, and an exchange reaction set as the product
      objective (model.c).

Output: ecTyrosol.mat in this directory. The base ecYeastGEM_batch.mat file is
never modified.
"""
import os
import numpy as np
import scipy.io as sio
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
# Override with environment variable ECYEASTGEM_BATCH if the base model lives elsewhere.
DEFAULT_BASE = os.path.expanduser(
    "~/Documents/CellFactory-ecYeastGEM/ModelFiles/ecYeastGEM_batch.mat"
)
BASE = os.environ.get("ECYEASTGEM_BATCH", DEFAULT_BASE)
OUT = os.path.join(HERE, "ecTyrosol.mat")

# Assumption A4: representative kcat for dedicated production arms (1/s)
KCAT_ARO10 = 1000.0
KCAT_ADH7 = 1000.0
S_ARO10 = 1.0 / (KCAT_ARO10 * 3600.0)
S_ADH7 = 1.0 / (KCAT_ADH7 * 3600.0)


def load(path):
    raw = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
    for k, v in raw.items():
        if not k.startswith("__") and hasattr(v, "_fieldnames"):
            return k, v
    raise RuntimeError("no struct in %s" % path)


def as_arr1d(x):
    if x is None:
        return np.array([], dtype=object)
    a = np.asarray(x)
    if a.ndim == 0:
        a = np.array([a.item()], dtype=object)
    return a.ravel()


def to_str(x):
    return "" if x is None else str(x)


def find_cyt_index(mets, metNames, metComps, compNames, target_name):
    comps = [to_str(c) for c in metComps]
    cytoplasm_id = None
    for i, n in enumerate(compNames):
        if to_str(n).lower() == "cytoplasm":
            cytoplasm_id = i + 1
            break
    for i, nm in enumerate(metNames):
        if to_str(nm) == target_name:
            c = comps[i] if i < len(comps) else ""
            try:
                if int(float(c)) == cytoplasm_id:
                    return i
            except Exception:
                if c.lower() in ("c", "cytoplasm"):
                    return i
    for i, nm in enumerate(metNames):
        if to_str(nm) == target_name:
            return i
    raise RuntimeError("metabolite %s not found" % target_name)


def find_met_by_id(mets, mid):
    for i, x in enumerate(mets):
        if to_str(x) == mid:
            return i
    raise RuntimeError("met id %s not found" % mid)


def cytoplasm_comp_token(metComps, compNames):
    for i, n in enumerate(compNames):
        if to_str(n).lower() == "cytoplasm":
            return i + 1
    return 1


def extracellular_comp_token(metComps, compNames):
    for i, n in enumerate(compNames):
        if to_str(n).lower() == "extracellular":
            return i + 1
    return 3


def main():
    _, m = load(BASE)
    print("Base model ecYeastGEM_batch: n_rxns=%d, n_mets=%d" % (m.S.shape[1], m.S.shape[0]))

    mets = [to_str(x) for x in as_arr1d(m.mets)]
    metNames = [to_str(x) for x in as_arr1d(m.metNames)]
    metFormulas = [to_str(x) for x in as_arr1d(m.metFormulas)] if hasattr(m, "metFormulas") else [""] * len(mets)
    metMiriams = list(as_arr1d(m.metMiriams)) if hasattr(m, "metMiriams") else [None] * len(mets)
    metComps = list(as_arr1d(m.metComps))
    compNames = [to_str(x) for x in as_arr1d(m.compNames)]
    cyt_token = cytoplasm_comp_token(metComps, compNames)
    ext_token = extracellular_comp_token(metComps, compNames)

    # Assumption A1/A2: 4-HPP is already in the batch model (tyrosine pathway)
    idx_HPP = find_cyt_index(mets, metNames, metComps, compNames, "3-(4-hydroxyphenyl)pyruvate")
    idx_H = find_met_by_id(mets, "s_0794")
    idx_CO2 = find_met_by_id(mets, "s_0456")
    idx_NADPH = find_met_by_id(mets, "s_1212")
    idx_NADP = find_met_by_id(mets, "s_1207")
    idx_pQ06408 = find_met_by_id(mets, "prot_Q06408")  # ARO10 (assumption A2)
    idx_pP25377 = find_met_by_id(mets, "prot_P25377")  # ADH7 (assumption A3)

    # Assumption A6: new tyrosol pool and export metabolites
    new_mets = [
        ("s_4hpaa_c", "4-hydroxyphenylacetaldehyde", cyt_token, "C8H8O2"),
        ("s_tyrosol_c", "tyrosol", cyt_token, "C8H10O2"),
        ("s_tyrosol_e", "tyrosol", ext_token, "C8H10O2"),
    ]
    new_idx = {}
    for mid, nm, ct, formula in new_mets:
        mets.append(mid)
        metNames.append(nm)
        metFormulas.append(formula)
        metComps.append(ct)
        metMiriams.append(None)
        new_idx[mid] = len(mets) - 1

    new_rxns = []

    # Assumption A2 + A4: ARO10-catalysed 4-HPP -> 4-HPAA
    new_rxns.append((
        "new_aro10_HPP",
        "4-hydroxyphenylpyruvate decarboxylase (ARO10)",
        "YDR380W", 0.0, 1000.0, "4.1.1.43",
        "sce00350  Tyrosine metabolism",
        {
            idx_H: -1.0,
            idx_HPP: -1.0,
            idx_pQ06408: -S_ARO10,
            new_idx["s_4hpaa_c"]: +1.0,
            idx_CO2: +1.0,
        },
    ))

    # Assumption A3 + A4: ADH7-catalysed 4-HPAA -> tyrosol
    new_rxns.append((
        "new_adh7_tyrosol",
        "4-hydroxyphenylacetaldehyde reductase (ADH7)",
        "YCR105W", 0.0, 1000.0, "1.1.1.90",
        "sce00350  Tyrosine metabolism",
        {
            idx_H: -1.0,
            idx_NADPH: -1.0,
            new_idx["s_4hpaa_c"]: -1.0,
            idx_pP25377: -S_ADH7,
            idx_NADP: +1.0,
            new_idx["s_tyrosol_c"]: +1.0,
        },
    ))

    # Assumption A6: transport and exchange
    new_rxns.append((
        "new_tyrosol_t", "tyrosol transport", "", 0.0, 1000.0, "",
        "sce04147  Exosome",
        {new_idx["s_tyrosol_c"]: -1.0, new_idx["s_tyrosol_e"]: +1.0},
    ))
    new_rxns.append((
        "new_tyrosol_ex", "tyrosol exchange", "", 0.0, 1000.0, "",
        "",
        {new_idx["s_tyrosol_e"]: -1.0},
    ))

    S = sp.csc_matrix(m.S) if not sp.issparse(m.S) else sp.csc_matrix(m.S)
    n_old_mets, n_old_rxns = S.shape
    S = sp.vstack([S, sp.csc_matrix((len(new_mets), n_old_rxns), dtype=S.dtype)]).tocsc()

    cols_data = []
    for *_, coldict in new_rxns:
        c = sp.lil_matrix((S.shape[0], 1), dtype=S.dtype)
        for r, v in coldict.items():
            c[r, 0] = v
        cols_data.append(c.tocsc())
    S = sp.hstack([S] + cols_data).tocsc()

    rxns = [to_str(x) for x in as_arr1d(m.rxns)]
    rxnNames = [to_str(x) for x in as_arr1d(m.rxnNames)]
    lb = list(np.asarray(m.lb).ravel())
    ub = list(np.asarray(m.ub).ravel())
    c_vec = [0.0] * n_old_rxns
    rev = list(np.asarray(m.rev).ravel()) if hasattr(m, "rev") else [0] * n_old_rxns
    grRules = [to_str(x) for x in as_arr1d(m.grRules)]
    subSystems = list(as_arr1d(m.subSystems)) if hasattr(m, "subSystems") else [""] * n_old_rxns
    eccodes = [to_str(x) for x in as_arr1d(m.eccodes)] if hasattr(m, "eccodes") else [""] * n_old_rxns
    rxnMiriams = list(as_arr1d(m.rxnMiriams)) if hasattr(m, "rxnMiriams") else [None] * n_old_rxns
    rxnReferences = [to_str(x) for x in as_arr1d(m.rxnReferences)] if hasattr(m, "rxnReferences") else [""] * n_old_rxns
    rxnConfidenceScores = list(np.asarray(m.rxnConfidenceScores).ravel()) if hasattr(m, "rxnConfidenceScores") else [0.0] * n_old_rxns

    objective_rxn_id = "new_tyrosol_ex"
    for rid, name, gr, _lb, _ub, ec, sub, _ in new_rxns:
        rxns.append(rid)
        rxnNames.append(name)
        lb.append(_lb)
        ub.append(_ub)
        rev.append(0)
        c_vec.append(1.0 if rid == objective_rxn_id else 0.0)
        grRules.append(gr)
        subSystems.append(sub)
        eccodes.append(ec)
        rxnMiriams.append(None)
        rxnReferences.append("")
        rxnConfidenceScores.append(2)

    genes = [to_str(x) for x in as_arr1d(m.genes)]
    rgm = sp.csc_matrix(m.rxnGeneMat) if not sp.issparse(m.rxnGeneMat) else sp.csc_matrix(m.rxnGeneMat)
    new_rgm_rows = sp.lil_matrix((len(new_rxns), len(genes)), dtype=rgm.dtype)
    for i, (_, _, gr, *_) in enumerate(new_rxns):
        if not gr:
            continue
        for tok in gr.replace("(", " ").replace(")", " ").replace(" or ", " ").replace(" and ", " ").split():
            if tok in genes:
                new_rgm_rows[i, genes.index(tok)] = 1
    rgm = sp.vstack([rgm, new_rgm_rows.tocsc()]).tocsc()

    b = list(np.asarray(m.b).ravel()) + [0.0] * len(new_mets)
    csense = list(as_arr1d(m.csense)) if hasattr(m, "csense") else None
    if csense is not None:
        csense = [to_str(x) for x in csense] + ["E"] * len(new_mets)

    raw = sio.loadmat(BASE, squeeze_me=False, struct_as_record=False)
    struct_var = next(k for k in raw if not k.startswith("__"))
    base_struct = raw[struct_var][0, 0]
    out = {fname: base_struct.__dict__[fname] for fname in base_struct._fieldnames}

    out["rxns"] = np.array(rxns, dtype=object).reshape(-1, 1)
    out["mets"] = np.array(mets, dtype=object).reshape(-1, 1)
    out["S"] = S
    out["lb"] = np.array(lb, dtype=float).reshape(-1, 1)
    out["ub"] = np.array(ub, dtype=float).reshape(-1, 1)
    out["rev"] = np.array(rev, dtype=float).reshape(-1, 1)
    out["c"] = np.array(c_vec, dtype=float).reshape(-1, 1)
    out["b"] = np.array(b, dtype=float).reshape(-1, 1)
    out["rxnNames"] = np.array(rxnNames, dtype=object).reshape(-1, 1)
    out["metNames"] = np.array(metNames, dtype=object).reshape(-1, 1)
    out["metFormulas"] = np.array(metFormulas, dtype=object).reshape(-1, 1)
    out["metComps"] = np.array(metComps, dtype=float).reshape(-1, 1)
    out["grRules"] = np.array(grRules, dtype=object).reshape(-1, 1)
    out["subSystems"] = np.array(subSystems, dtype=object).reshape(-1, 1)
    out["eccodes"] = np.array(eccodes, dtype=object).reshape(-1, 1)
    out["rxnConfidenceScores"] = np.array(rxnConfidenceScores, dtype=float).reshape(-1, 1)
    out["rxnReferences"] = np.array(rxnReferences, dtype=object).reshape(-1, 1)
    out["rxnGeneMat"] = rgm
    out["description"] = (
        "ecTyrosol: ecYeastGEM_batch extended for tyrosol production via the "
        "Ehrlich pathway (ARO10, ADH7) with enzyme-constrained arms (kcat=1000/s)."
    )

    empty = np.empty((0, 0), dtype=object)
    if "metCharges" in out:
        old = np.asarray(out["metCharges"]).reshape(-1, 1).astype(float)
        out["metCharges"] = np.vstack([old, np.array([[np.nan]] * len(new_mets))])
    if "concs" in out:
        old = np.asarray(out["concs"]).reshape(-1, 1).astype(float)
        out["concs"] = np.vstack([old, np.array([[np.nan]] * len(new_mets))])
    if "metMiriams" in out:
        old_arr = np.asarray(out["metMiriams"], dtype=object).reshape(-1, 1)
        new_rows = np.array([[empty]] * len(new_mets), dtype=object).reshape(-1, 1)
        out["metMiriams"] = np.vstack([old_arr, new_rows])
    if "rxnMiriams" in out:
        old_arr = np.asarray(out["rxnMiriams"], dtype=object).reshape(-1, 1)
        new_rows = np.array([[empty]] * len(new_rxns), dtype=object).reshape(-1, 1)
        out["rxnMiriams"] = np.vstack([old_arr, new_rows])
    if csense is not None:
        out["csense"] = np.array(csense, dtype=object).reshape(-1, 1)

    sio.savemat(OUT, {struct_var: out}, do_compression=False)
    print("Saved %s (%d rxns, %d mets)" % (OUT, len(rxns), len(mets)))
    print("Product objective: %s" % objective_rxn_id)


if __name__ == "__main__":
    main()
