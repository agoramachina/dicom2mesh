# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pydicom>=3",
#     "python-gdcm",
#     "numpy",
#     "scipy",
#     "scikit-image",
#     "trimesh",
#     "pymeshlab",
#     "nibabel",
#     "brainchop",
# ]
# ///
"""
brain2mesh - turn a T1-weighted brain MRI (DICOM) into a 3D brain surface.

    uv run brain2mesh.py SCANS_FOLDER                      # list the series found
    uv run brain2mesh.py SCANS_FOLDER -s 10 -o brain.stl   # convert series #10

Uses brainchop's MeshNet models to strip the skull and label brain tissue,
then places the surface using the T1 intensities themselves. Works best on a
3D T1 (e.g. MPRAGE). Proof-of-concept; not for clinical use.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import ndimage

import dicom2mesh as d

# brainchop's "subcortical" model uses SynthSeg-style labels with left and
# right merged. The model ships without names; these two were identified by
# size, position and intensity (cortex ribbon, dark ventricles).
CORTEX, LATERAL_VENTRICLE = 2, 3
LPS_TO_RAS = np.diag([-1.0, -1.0, 1.0])


def save_nifti(vol, A, origin, path: Path):
    """Write the volume with its voxel->patient transform (NIfTI uses RAS, DICOM LPS)."""
    aff = np.eye(4)
    aff[:3, :3] = LPS_TO_RAS @ A
    aff[:3, 3] = LPS_TO_RAS @ origin
    img = nib.Nifti1Image(vol.astype(np.float32), aff)
    img.set_qform(aff, 1)
    img.set_sform(aff, 1)
    nib.save(img, path)


def brainchop(args, output: Path, device):
    """Run the brainchop CLI from this environment and fail loudly."""
    here = Path(sys.executable).parent
    exe = shutil.which("brainchop", path=str(here)) or shutil.which("brainchop")
    if exe is None:
        sys.exit("brainchop executable not found (is it installed in this environment?)")
    env = dict(os.environ, **({"DEV": device} if device else {}))
    r = subprocess.run([exe, *args, "-o", str(output), "--no-optimize"], capture_output=True,
                       encoding="utf-8", errors="replace", env=env)
    # Judge by the output file: brainchop's exit status isn't reliable on errors.
    if r.returncode != 0 or not output.exists():
        lines = (r.stderr + r.stdout).strip().splitlines()
        err = next((ln for ln in reversed(lines) if "Error" in ln), lines[-1] if lines else "(no output)")
        hint = ("\n  out of GPU memory: close other GPU apps, or rerun with --device CPU (slow)"
                if "MemoryError" in err else "")
        sys.exit(f"brainchop failed: {err}{hint}")


def label_brain(vol, A, origin, work: Path, device):
    """Skull-strip, then label tissue. Returns labels on the original voxel grid."""
    t1 = work / "t1.nii.gz"
    save_nifti(vol, A, origin, t1)
    # Stripping first lets --crop shrink the second pass to the brain, which
    # roughly halves GPU memory for the (larger) labelling model.
    d.log("brainchop: skull stripping (mindgrab)")
    brainchop([str(t1), "-m", "mindgrab", "--crop", "--inverse-conform"],
              work / "stripped.nii.gz", device)
    d.log("brainchop: tissue labels (subcortical)")
    brainchop([str(work / "stripped.nii.gz"), "-m", "subcortical", "--crop", "--inverse-conform"],
              work / "labels.nii.gz", device)
    labels = nib.load(work / "labels.nii.gz")
    if labels.shape != vol.shape:
        sys.exit(f"brainchop returned a {labels.shape} grid, expected {vol.shape}")
    return np.asarray(labels.dataobj).astype(np.uint8)


def pial_field(t1, labels, spacing, sigma_mm, reach):
    """Labels decide *what* is brain; T1 intensity decides *where* its surface is.

    Contouring binary labels gives a lumpy, voxel-noisy cortex. Contouring
    the T1 itself (between grey matter and CSF brightness) lets marching
    cubes place the pial surface at sub-voxel positions - but only near the
    labelled boundary, so enhancing vessels and dura can't join the brain.
    """
    cc, _ = ndimage.label(labels > 0)
    sizes = np.bincount(cc.ravel())
    sizes[0] = 0
    brain = ndimage.binary_fill_holes(cc == sizes.argmax())

    t1 = ndimage.gaussian_filter(t1, sigma_mm / np.asarray(spacing))
    ventricles = ndimage.binary_erosion(labels == LATERAL_VENTRICLE)
    if not (labels == CORTEX).any() or not ventricles.any():
        sys.exit("segmentation found no cortex or ventricles; is this a brain MRI?")
    gm = float(np.median(t1[labels == CORTEX]))
    # Eroded ventricle interior = pure CSF, free of partial-volume edges and
    # the (possibly contrast-enhanced) choroid plexus rim.
    csf = float(np.median(t1[ventricles]))
    if not csf < gm:
        sys.exit(f"CSF ({csf:.0f}) is not darker than grey matter ({gm:.0f}); "
                 "this needs a T1-weighted series (e.g. MPRAGE)")
    level = (gm + csf) / 2
    d.log(f"pial level: grey matter {gm:.0f}, CSF {csf:.0f} -> {level:.0f}")

    core = ndimage.binary_erosion(brain, iterations=2)           # always inside
    # (scipy treats iterations=0 as "dilate until nothing changes", so guard it)
    grown = ndimage.binary_dilation(brain, iterations=reach) if reach > 0 else brain
    band = grown & ~core
    field = np.zeros_like(t1)
    field[core] = 2 * gm
    field[band] = t1[band]
    field[band & ~brain & (t1 > gm)] = 0                         # vessels, dura
    return field, level


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=Path,
                   help="folder with DICOM files (a copied disc, or one series' folder)")
    p.add_argument("-s", "--series", type=int,
                   help="which series to convert, by its number in the listing")
    p.add_argument("-o", "--out", type=Path, action="append",
                   help=f"output file, repeatable; format from extension ({' '.join(d.FORMATS)})")
    p.add_argument("--sigma", type=float, default=0.7,
                   help="gaussian pre-smoothing of the T1 in mm (default 0.7)")
    p.add_argument("--smooth", type=int, default=15, help="taubin iterations (default 15)")
    p.add_argument("--faces", type=int, default=600_000,
                   help="decimate to about N faces, 0 = off (default 600000)")
    p.add_argument("--reach", type=int, default=1,
                   help="voxels the surface may move outward from the labels (default 1)")
    p.add_argument("--device", help="force brainchop's compute device, e.g. CPU (default: auto)")
    p.add_argument("--keep", type=Path,
                   help="keep brainchop's intermediate NIfTI files in this folder")
    a = p.parse_args()

    if not a.input.is_dir():
        sys.exit(f"not a folder: {a.input}")
    for out in a.out or []:
        if out.suffix.lower() not in d.FORMATS:
            sys.exit(f"unsupported output format '{out.suffix}'; use one of: {' '.join(d.FORMATS)}")

    series = d.find_series(a.input)
    if not series:
        sys.exit(f"no DICOM image series found under {a.input}")
    if not a.out:
        d.print_series(series, a.input)
        print("Pick a 3D T1 series (e.g. MPRAGE), then:  "
              "uv run brain2mesh.py <folder> -s <#> -o brain.stl")
        return
    if a.series is None and len(series) > 1:
        d.print_series(series, a.input)
        sys.exit("More than one series here: add -s <#> to pick one.")
    if a.series is not None and not 1 <= a.series <= len(series):
        sys.exit(f"--series must be between 1 and {len(series)}")

    ds, paths = series[(a.series or 1) - 1]
    if str(ds.get("Modality", "")) != "MR":
        sys.exit(f"series is {ds.get('Modality', '?')}, not MR; brain2mesh needs a T1 MRI")
    vol, A, origin = d.load_series(paths)
    spacing = np.linalg.norm(A, axis=0)

    if a.keep:
        a.keep.mkdir(parents=True, exist_ok=True)
        labels = label_brain(vol, A, origin, a.keep, a.device)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            labels = label_brain(vol, A, origin, Path(tmp), a.device)

    field, level = pial_field(vol, labels, spacing, a.sigma, a.reach)
    del vol, labels
    # Same cleanup as CT: keep one big piece, fill every enclosed pocket.
    field, bg = d.segment(field, spacing, level, 0, 50_000, 1e9, None)
    mesh = d.extract_surface(field, level, bg, A, origin)
    del field
    mesh = d.cleanup(mesh, a.smooth, a.faces)

    d.log(f"final: {len(mesh.faces):,} faces, watertight={mesh.is_watertight}, "
          f"volume={mesh.volume / 1000:.0f} cm^3, "
          f"extent={np.ptp(mesh.vertices, axis=0).round(1)} mm (L-R, A-P, S-I)")
    for out in a.out:
        d.export(mesh, out)


if __name__ == "__main__":
    main()
