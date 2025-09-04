# ======== PARAMETERS ========
Nmol = 750              # number of H2 molecules (atoms = 2*Nmol)
bond_ang = 0.74         # H-H bond length
spacing_ang = 2.5       # spacing between molecules
basis_name = "6-31G"
xc_functional = "lda"
charge = 0
spin_multiplicity = 1
threshold = 0        # sparsification threshold
save_prefix = f"h2_chain_{Nmol}mol_{basis_name}_{xc_functional}"
# ============================

import numpy as np
import time

try:
    from pyscf import gto, dft
except ModuleNotFoundError as e:
    raise SystemExit("PySCF missing. Install with: pip install pyscf") from e

try:
    from scipy.io import savemat
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False


def build_h2_chain(Nmol, bond_ang, spacing_ang):
    atoms = []
    half = bond_ang / 2.0
    for i in range(Nmol):
        zc = i * spacing_ang
        atoms.append(("H", (0.0, 0.0, zc - half)))
        atoms.append(("H", (0.0, 0.0, zc + half)))
    return atoms


def run_rks(atoms):
    mol = gto.Mole()
    mol.atom = atoms
    mol.unit = "Angstrom"
    mol.basis = basis_name
    mol.charge = charge
    mol.spin = 0 if spin_multiplicity == 1 else (spin_multiplicity - 1)
    mol.build()

    mf = dft.RKS(mol)
    mf.xc = xc_functional
    mf.conv_tol = 1e-8
    mf.max_cycle = 150
    mf.init_guess = 'atom'
    mf.damp = 0.2
    mf.diis_start_cycle = 5

    e_tot = mf.kernel()

    H_core = mf.get_hcore()
    S = mf.get_ovlp()
    D = mf.make_rdm1()
    V_eff = mf.get_veff(mf.mol, D)
    H_dft = H_core + V_eff
    return mol, mf, e_tot, H_dft, S


def homolumo_gap(mf):
    mo_e = mf.mo_energy
    nocc = mf.mol.nelectron // 2
    if nocc >= mo_e.size:
        return mo_e[nocc - 1], None, None
    HOMO = mo_e[nocc - 1]
    LUMO = mo_e[nocc]
    return HOMO, LUMO, LUMO - HOMO


def relative_spectral_gap(mf):
    mo_e = mf.mo_energy
    nocc = mf.mol.nelectron // 2
    N = mo_e.size
    if nocc >= N:
        return None
    gap = mo_e[nocc] - mo_e[nocc - 1]
    width = mo_e[N - 1] - mo_e[0]
    return gap / width if width != 0 else None


def sparsify_and_stats(M, thr):
    A = M.copy()
    A[np.abs(A) < thr] = 0.0
    nnz = np.count_nonzero(A)
    frac = nnz / A.size
    return A, nnz, frac


def main():
    atoms = build_h2_chain(Nmol, bond_ang, spacing_ang)
    Nat = len(atoms)
    Nelec = Nat  # 1e per H atom
    print(f"Built H2 chain: Nmol = {Nmol} (atoms = {Nat}, electrons = {Nelec})")

    mol, mf, e_tot, H_dft, S = run_rks(atoms)
    nbf = H_dft.shape[0]
    print(f"AO basis size nbf = {nbf}")

    HOMO, LUMO, gap = homolumo_gap(mf)
    print(f"Total energy: {e_tot:.8f} Ha")
    if LUMO is None:
        print(f"HOMO = {HOMO:.6f} Ha, LUMO = None (no virtuals), Gap = undefined")
    else:
        print(f"HOMO = {HOMO:.6f} Ha, LUMO = {LUMO:.6f} Ha, "
              f"Gap = {gap:.6f} Ha ({gap*27.2114:.2f} eV)")

    # Relative spectral gap
    mf.kernel()
    rel_gap = relative_spectral_gap(mf)
    if rel_gap is not None:
        print(f"Relative spectral gap: {rel_gap:.6f}")
    else:
        print("Relative spectral gap: undefined")

    # Save eigenvalues
    mo_e = mf.mo_energy
    np.save(f"{save_prefix}_eigs.npy", mo_e)
    np.savetxt(f"{save_prefix}_eigs.csv", mo_e, delimiter=",")
    print(f"Saved eigenvalues to {save_prefix}_eigs.npy/.csv")

    # Sparsify
    H_thr, nnzH, fracH = sparsify_and_stats(H_dft, threshold)
    S_thr, nnzS, fracS = sparsify_and_stats(S, threshold)
    print(f"Sparsity after |x|<{threshold:g} --> 0:")
    print(f"  H_DFT: nnz={nnzH:,}  density={100*fracH:.2f}%")
    print(f"  S    : nnz={nnzS:,}  density={100*fracS:.2f}%")

    # Save matrices
    np.save(f"{save_prefix}_H_DFT.npy", H_dft)
    np.save(f"{save_prefix}_S.npy", S)
    np.savetxt(f"{save_prefix}_H_DFT.csv", H_dft, delimiter=",")
    np.savetxt(f"{save_prefix}_S.csv", S, delimiter=",")
    np.save(f"{save_prefix}_H_DFT_thr{threshold}.npy", H_thr)
    np.save(f"{save_prefix}_S_thr{threshold}.npy", S_thr)

    if HAVE_SCIPY:
        savemat(f"{save_prefix}_matrices.mat",
                {"H_DFT": H_dft, "S": S,
                 "H_DFT_thr": H_thr, "S_thr": S_thr,
                 "eigs_Ha": mo_e,
                 "HOMO": HOMO, "LUMO": LUMO,
                 "gap_Ha": (gap if gap is not None else np.nan),
                 "gap_eV": (gap*27.2114 if gap is not None else np.nan),
                 "rel_gap": (rel_gap if rel_gap is not None else np.nan)})
        print(f"Wrote {save_prefix}_matrices.mat")


if __name__ == "__main__":
    start_time = time.time()
    main()
    print("Program execution time: %s seconds." % (time.time() - start_time))
