#!/usr/bin/env python3
"""
Read g-tensor or hyperfine A-tensors directly from a Quantum ESPRESSO
GIPAW (gipaw.x) output, diagonalize them, and print principal axes plus
VESTA-ready directions.

For g-tensor:
    - reads the LAST "Delta_g total ..." 3x3 tensor in the output
      (there may be both "SOO a la Paratec" and "SOO as in Eq.(7)";
       the last one is selected).

For hyperfine:
    - reads the LAST "total dipolar (symmetrized)" 3x3 tensor for each atom
    - reads "Fermi contact in <unit>" and takes the "total" column
    - constructs the full hyperfine tensor:
          A_total = A_dipolar(symmetrized) + A_FC(total) * I
    - by default processes atoms 1..7.

Notes on lattice/VESTA:
    gipaw.x output itself does not print the lattice vectors in the supplied
    examples. Therefore this script defaults to Cartesian unit vectors:
        a=(1,0,0), b=(0,1,0), c=(0,0,1)
    This is sufficient when the simulation-cell axes are aligned with Cartesian
    x,y,z (as for a conventional cubic diamond supercell).

    Preferred for a general cell: supply the original pw.x input:
        --qe-input scf.in

    The script reads CELL_PARAMETERS automatically. Manual lattice vectors
    remain available:
        --a ax ay az --b bx by bz --c cx cy cz

    Crystallographic z/y definitions can be supplied as integer combinations
    of a,b,c:
        --iz ia ib ic
        --iy ja jb jc
    Defaults are z=c and y=b.

Examples:
    python diag_tensor_gipaw.py g-tensor.out
    python diag_tensor_gipaw.py hyperfine.out
    python diag_tensor_gipaw.py hyperfine.out --atoms 1-7
    python diag_tensor_gipaw.py hyperfine.out --qe-input scf.in

    python diag_tensor_gipaw.py g-tensor.out \
        --a 10 0 0 --b 0 10 0 --c 0 0 10

Dependencies:
    numpy
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
import numpy as np

PI = math.pi
FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def fnum(s: str) -> float:
    return float(s.replace("D", "E").replace("d", "e"))


def detect_job(text: str) -> str:
    m = re.search(r"job\s*=\s*['\"]([^'\"]+)['\"]", text, re.I)
    if m:
        job = m.group(1).strip().lower()
        if "g_tensor" in job:
            return "g"
        if "hyperfine" in job:
            return "A"
    if "Delta_g total" in text:
        return "g"
    if "Fermi contact" in text and "dipolar" in text:
        return "A"
    raise ValueError("Could not determine whether this is a g_tensor or hyperfine GIPAW output.")


def parse_last_g_tensor(text: str):
    """Return the last Delta_g total 3x3 tensor and its heading."""
    lines = text.splitlines()
    hits = [i for i, line in enumerate(lines) if "Delta_g total" in line]
    if not hits:
        raise ValueError('No "Delta_g total" block found.')

    i = hits[-1]
    rows = []
    j = i + 1
    while j < len(lines) and len(rows) < 3:
        vals = re.findall(FLOAT_RE, lines[j])
        if len(vals) >= 3:
            rows.append([fnum(x) for x in vals[-3:]])
        j += 1
    if len(rows) != 3:
        raise ValueError("Could not read 3 rows of the last Delta_g total tensor.")
    return np.array(rows, float), lines[i].strip()


def _parse_tensor_block(lines, heading_index, max_atoms=None):
    """Parse atom-labelled 3x3 tensors after a heading."""
    tensors = {}
    labels = {}
    i = heading_index + 1

    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.startswith("-----") and i > heading_index + 1:
            break
        if "Initialization:" in line or "Fermi contact" in line:
            break

        # Expected row:
        # Sn  1   x y z
        # C1  2   x y z
        parts = line.split()
        if len(parts) >= 5:
            try:
                atom_index = int(parts[1])
                row = [fnum(parts[2]), fnum(parts[3]), fnum(parts[4])]
            except (ValueError, IndexError):
                i += 1
                continue

            label = parts[0]
            rows = [row]
            k = i + 1
            while k < len(lines) and len(rows) < 3:
                p = lines[k].split()
                if len(p) >= 5:
                    try:
                        idx2 = int(p[1])
                        if idx2 != atom_index:
                            break
                        rows.append([fnum(p[2]), fnum(p[3]), fnum(p[4])])
                    except (ValueError, IndexError):
                        break
                else:
                    if lines[k].strip():
                        break
                k += 1

            if len(rows) == 3:
                tensors[atom_index] = np.array(rows, float)
                labels[atom_index] = label
                if max_atoms is not None and len(tensors) >= max_atoms:
                    pass
                i = k
                continue
        i += 1

    return tensors, labels


def parse_last_symmetrized_dipolar(text: str):
    lines = text.splitlines()
    hits = [i for i, line in enumerate(lines)
            if "total dipolar (symmetrized)" in line.lower()]
    if not hits:
        raise ValueError('No "total dipolar (symmetrized)" block found.')
    return _parse_tensor_block(lines, hits[-1])


def parse_fermi_contact_total(text: str):
    """
    Parse the final Fermi-contact table.
    Returns {atom_index: total_FC}, labels, unit.
    """
    lines = text.splitlines()
    hits = [i for i, line in enumerate(lines) if "Fermi contact in" in line]
    if not hits:
        raise ValueError('No "Fermi contact in ..." table found.')
    h = hits[-1]
    heading = lines[h]
    m = re.search(r"Fermi contact in\s+(\S+)", heading, re.I)
    unit = m.group(1) if m else "unknown"

    vals = {}
    labels = {}
    i = h + 1
    # skip column heading(s)
    while i < len(lines):
        p = lines[i].split()
        if len(p) >= 6:
            try:
                idx = int(p[1])
                total = fnum(p[-1])
                vals[idx] = total
                labels[idx] = p[0]
            except ValueError:
                pass
        elif vals and (not lines[i].strip() or "Initialization:" in lines[i]):
            if "Initialization:" in lines[i]:
                break
        if "Initialization:" in lines[i]:
            break
        i += 1

    if not vals:
        raise ValueError("Could not parse Fermi-contact total column.")
    return vals, labels, unit


def parse_atoms(spec: str) -> list[int]:
    result = []
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            a, b = item.split("-", 1)
            a, b = int(a), int(b)
            step = 1 if b >= a else -1
            result.extend(range(a, b + step, step))
        else:
            result.append(int(item))
    # preserve order, remove duplicates
    return list(dict.fromkeys(result))


BOHR_TO_ANG = 0.529177210903


def _qe_scalar(text: str, name: str):
    """Read a scalar QE namelist parameter such as A=..., celldm(1)=... ."""
    pat = rf"(?im)^\s*{re.escape(name)}\s*=\s*({FLOAT_RE})"
    m = re.search(pat, text)
    return fnum(m.group(1)) if m else None


def read_qe_cell(path: Path):
    """
    Read lattice vectors from a Quantum ESPRESSO pw.x input.

    Supported:
      1) CELL_PARAMETERS angstrom
      2) CELL_PARAMETERS bohr
      3) CELL_PARAMETERS alat
         with alat taken from A (angstrom) or celldm(1) (bohr)
      4) CELL_PARAMETERS with no explicit unit:
         QE convention is alat; A/celldm(1) is therefore required here.

    Returns:
        aa, ab, ac in Angstrom, and a human-readable description.

    For ibrav != 0 without CELL_PARAMETERS, this function currently raises
    an informative error rather than silently reconstructing a Bravais cell.
    """
    text = path.read_text(errors="replace")
    lines = text.splitlines()

    # Find the last CELL_PARAMETERS block, matching QE behavior if a file
    # accidentally contains more than one such block.
    hits = [i for i, line in enumerate(lines)
            if re.match(r"^\s*CELL_PARAMETERS\b", line, re.I)]

    A_ang = _qe_scalar(text, "A")
    celldm1_bohr = _qe_scalar(text, "celldm(1)")
    ibrav_val = _qe_scalar(text, "ibrav")
    ibrav = int(round(ibrav_val)) if ibrav_val is not None else None

    if hits:
        i = hits[-1]
        header = lines[i].strip()

        unit_m = re.search(
            r"CELL_PARAMETERS\s*(?:\(\s*)?(angstrom|bohr|alat)(?:\s*\))?",
            header, re.I
        )
        unit = unit_m.group(1).lower() if unit_m else "alat"

        rows = []
        for j in range(i + 1, min(i + 4, len(lines))):
            vals = re.findall(FLOAT_RE, lines[j])
            if len(vals) < 3:
                raise ValueError(
                    f"Could not read three lattice-vector components from "
                    f"{path}, line {j+1}: {lines[j]!r}"
                )
            rows.append([fnum(x) for x in vals[:3]])

        cell = np.array(rows, dtype=float)

        if unit == "angstrom":
            factor = 1.0
            source = "CELL_PARAMETERS angstrom"
        elif unit == "bohr":
            factor = BOHR_TO_ANG
            source = "CELL_PARAMETERS bohr"
        else:  # alat
            if A_ang is not None:
                factor = A_ang
                source = f"CELL_PARAMETERS alat, A={A_ang:g} angstrom"
            elif celldm1_bohr is not None:
                factor = celldm1_bohr * BOHR_TO_ANG
                source = (
                    "CELL_PARAMETERS alat, "
                    f"celldm(1)={celldm1_bohr:g} bohr"
                )
            else:
                raise ValueError(
                    "CELL_PARAMETERS is in alat units, but neither A nor "
                    "celldm(1) was found in the QE input, so alat cannot be "
                    "converted to an absolute lattice."
                )

        cell *= factor
        return cell[0], cell[1], cell[2], source

    # No CELL_PARAMETERS.
    if ibrav == 0:
        raise ValueError(
            f"{path}: ibrav=0 but no CELL_PARAMETERS block was found."
        )

    raise ValueError(
        f"{path}: no CELL_PARAMETERS block found"
        + (f" (ibrav={ibrav})" if ibrav is not None else "")
        + ". Automatic reconstruction for ibrav != 0 is not implemented; "
          "either provide a QE input containing CELL_PARAMETERS or use "
          "--a/--b/--c explicitly."
    )


def norm(v):
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("Cannot normalize zero vector.")
    return v / n


def symmetrize(A):
    return 0.5 * (A + A.T)


def crystallographic_axes(aa, ab, ac, iz, iy):
    az = iz[0]*aa + iz[1]*ab + iz[2]*ac
    ay = iy[0]*aa + iy[1]*ab + iy[2]*ac
    ax = np.cross(az, ay)
    ax, ay, az = norm(ax), norm(ay), norm(az)
    return ax, ay, az


def theta_phi(q, az, ay, ax):
    costheta = float(np.clip(np.dot(q, az), -1.0, 1.0))
    theta = math.acos(costheta)
    qxy = q - costheta*az
    qn = np.linalg.norm(qxy)
    if qn < 1e-14:
        phi = 0.0
    else:
        cosphi = float(np.clip(np.dot(qxy, ax)/qn, -1.0, 1.0))
        phi = math.acos(cosphi)
        if np.dot(qxy, ay) < 0:
            phi = -phi
    return theta, phi


def reduced_angles(theta, phi):
    if theta > 0.5*PI:
        theta = PI - theta
    if phi > 0.5*PI:
        phi = -(PI - phi)
    elif phi < -0.5*PI:
        phi = PI + phi
    return theta, phi


def lattice_components(q, aa, ab, ac):
    basis = np.column_stack((aa, ab, ac))
    return np.linalg.solve(basis, q)


def print_matrix(A, fmt="12.6f"):
    for row in A:
        print(" ".join(f"{x:{fmt}}" for x in row))


def analyze_tensor(A, name, tensor_type, aa, ab, ac, iz, iy, unit=""):
    Araw = np.array(A, float)
    A = symmetrize(Araw)

    ax, ay, az = crystallographic_axes(aa, ab, ac, iz, iy)
    W, Q = np.linalg.eigh(A)

    suffix = f" {unit}" if unit else ""

    print("\n" + "="*88)
    print(name)
    print("="*88)
    print("\nTensor used for diagonalization:")
    print_matrix(A)

    print("\nPrincipal values:" + "".join(f"{x:14.6f}" for x in W) + suffix)

    print("\nPrincipal axes in Cartesian coordinates:")
    for j in range(3):
        q = Q[:, j]
        print(f"  {j+1}: {q[0]:14.8f} {q[1]:14.8f} {q[2]:14.8f}"
              f"   |q|={np.linalg.norm(q):.12f}")

    print("\nPrincipal-axis angles and VESTA vectors:")
    for j in range(3):
        q = Q[:, j]
        th, ph = theta_phi(q, az, ay, ax)
        thr, phr = reduced_angles(th, ph)
        uvw = lattice_components(q, aa, ab, ac)

        # Keep the historical Fortran visualization convention:
        # g: scale by eigenvalue/10000; A: use fixed length 1.1.
        if tensor_type == "g":
            vesta = uvw * (W[j] / 10000.0)
        else:
            vesta = uvw * 1.1

        print(f"\n  axis {j+1}, principal value = {W[j]:.6f}{suffix}")
        print(f"    theta = {math.degrees(th):8.3f} deg"
              f" ({math.degrees(thr):8.3f})")
        print(f"    phi   = {math.degrees(ph):8.3f} deg"
              f" ({math.degrees(phr):8.3f})")
        print(f"    uvw   = {uvw[0]:14.8f} {uvw[1]:14.8f} {uvw[2]:14.8f}")
        print("    VESTA:")
        print(f"      { vesta[0]:14.7f} { vesta[1]:14.7f} { vesta[2]:14.7f}")
        print(f"      {-vesta[0]:14.7f} {-vesta[1]:14.7f} {-vesta[2]:14.7f}")

    print("\nOrthogonality checks:")
    print(f"  q1.q2 = {np.dot(Q[:,0],Q[:,1]): .6e}")
    print(f"  q2.q3 = {np.dot(Q[:,1],Q[:,2]): .6e}")
    print(f"  q3.q1 = {np.dot(Q[:,2],Q[:,0]): .6e}")

    return W, Q


def run(args):
    text = args.input.read_text(errors="replace")
    job = args.type or detect_job(text)

    if args.qe_input is not None:
        aa, ab, ac, cell_source = read_qe_cell(args.qe_input)
    else:
        aa = np.array(args.a, float)
        ab = np.array(args.b, float)
        ac = np.array(args.c, float)
        cell_source = "command line/default Cartesian cell"

    iz = np.array(args.iz, int)
    iy = np.array(args.iy, int)

    print(f"Input file: {args.input}")
    print(f"Detected tensor type: {'g-tensor' if job == 'g' else 'hyperfine A-tensor'}")
    if args.qe_input is not None:
        print(f"QE input: {args.qe_input}")
    print(f"Cell source: {cell_source}")
    print("\nLattice vectors used for VESTA conversion (Angstrom):")
    print(f"  a = {aa}")
    print(f"  b = {ab}")
    print(f"  c = {ac}")
    print(f"  Iz = {iz}   (crystallographic z)")
    print(f"  Iy = {iy}   (crystallographic y)")

    if job == "g":
        A, heading = parse_last_g_tensor(text)
        print(f"\nSelected GIPAW block: {heading}")
        analyze_tensor(
            A, "EPR Delta_g total tensor", "g",
            aa, ab, ac, iz, iy, unit="ppm"
        )
        return

    dip, dip_labels = parse_last_symmetrized_dipolar(text)
    fc, fc_labels, unit = parse_fermi_contact_total(text)

    atoms = parse_atoms(args.atoms)
    print("\nHyperfine construction:")
    print("  A_total = A_dipolar(total, symmetrized) + A_Fermi-contact(total) * I")
    print(f"  Requested atoms: {atoms}")
    print(f"  Output unit: {unit}")

    for idx in atoms:
        if idx not in dip:
            print(f"\nWARNING: atom {idx}: dipolar tensor not found; skipped.")
            continue
        if idx not in fc:
            print(f"\nWARNING: atom {idx}: Fermi-contact total not found; skipped.")
            continue

        label = fc_labels.get(idx, dip_labels.get(idx, "?"))
        A_dip = dip[idx]
        A_fc = fc[idx]
        A_total = A_dip + np.eye(3)*A_fc

        print("\n" + "#"*88)
        print(f"Atom {idx}: {label}")
        print("#"*88)
        print(f"\nFermi contact (total) = {A_fc:.6f} {unit}")
        print("\nDipolar tensor (total, symmetrized):")
        print_matrix(A_dip)
        print("\nFull hyperfine tensor A_total:")
        print_matrix(A_total)

        analyze_tensor(
            A_total,
            f"Hyperfine A-tensor: atom {idx} ({label})",
            "A", aa, ab, ac, iz, iy, unit=unit
        )


def parse_args():
    p = argparse.ArgumentParser(
        description="Read GIPAW g/hyperfine tensors directly from gipaw.x output."
    )
    p.add_argument("input", type=Path, help="gipaw.x output file")
    p.add_argument(
        "--type", choices=("g", "A"), default=None,
        help="force tensor type; normally detected automatically"
    )
    p.add_argument(
        "--atoms", default="1-7",
        help="hyperfine atom indices, e.g. 1-7 or 1,2,5 (default: 1-7)"
    )
    p.add_argument(
        "--qe-input", type=Path, default=None,
        help="pw.x SCF input; read CELL_PARAMETERS automatically for VESTA conversion"
    )
    p.add_argument("--a", nargs=3, type=float, default=(1.0,0.0,0.0),
                   metavar=("AX","AY","AZ"),
                   help="lattice vector a; default Cartesian x")
    p.add_argument("--b", nargs=3, type=float, default=(0.0,1.0,0.0),
                   metavar=("BX","BY","BZ"),
                   help="lattice vector b; default Cartesian y")
    p.add_argument("--c", nargs=3, type=float, default=(0.0,0.0,1.0),
                   metavar=("CX","CY","CZ"),
                   help="lattice vector c; default Cartesian z")
    p.add_argument("--iz", nargs=3, type=int, default=(0,0,1),
                   metavar=("IA","IB","IC"),
                   help="crystallographic z as integer combination of a,b,c")
    p.add_argument("--iy", nargs=3, type=int, default=(0,1,0),
                   metavar=("JA","JB","JC"),
                   help="crystallographic y as integer combination of a,b,c")
    return p.parse_args()


if __name__ == "__main__":
    run(parse_args())
