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
# ]
# ///
"""
dicom2mesh - turn a DICOM series (e.g. a head CT) into a 3D mesh.

    uv run dicom2mesh.py SCANS_FOLDER                       # list the series found
    uv run dicom2mesh.py SCANS_FOLDER -s 3 -o skull.stl     # convert series #3

Proof-of-concept. Threshold-based, so it works well for CT bone (Hounsfield
units are calibrated), and poorly for MRI (intensities are not). Not for
clinical use.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pydicom
import pymeshlab
import trimesh
from scipy import ndimage
from skimage import measure

T0 = time.time()


def log(msg):
    # ASCII only: Windows consoles and redirected output may not be UTF-8.
    print(f"[{time.time() - T0:6.1f}s] {msg}", file=sys.stderr)


# ── 1. find & load ────────────────────────────────────────────────────────────

def find_series(root: Path):
    """Group every DICOM image under `root` by series, reading headers only."""
    files, first, empty = {}, {}, 0
    for f in sorted(p for p in root.rglob("*") if p.is_file()):
        if f.stat().st_size == 0:
            empty += 1
            continue
        try:
            ds = pydicom.dcmread(f, stop_before_pixels=True)
        except Exception:
            continue                              # viewer .exe, autorun.inf, ...
        if "ImagePositionPatient" in ds and "Rows" in ds:
            files.setdefault(ds.SeriesInstanceUID, []).append(f)
            first.setdefault(ds.SeriesInstanceUID, ds)
    if empty:
        log(f"warning: skipped {empty} empty (0-byte) files -- incomplete copy from the disc?")
    series = [(first[uid], fs) for uid, fs in files.items()]
    series.sort(key=lambda s: (str(s[1][0].parent), int(s[0].get("SeriesNumber") or 0)))
    return series


def print_series(series, root: Path):
    print(f"\nFound {len(series)} image series under {root}:\n")
    print(f"  {'#':>2}  {'slices':>6}  mod  {'slice':>6}  description / folder")
    for i, (ds, fs) in enumerate(series, 1):
        thick = float(ds.get("SliceThickness") or 0)
        hint = "   <- good for 3D" if len(fs) >= 50 and 0 < thick <= 1.5 else ""
        print(f"  {i:>2}  {len(fs):>6}  {str(ds.get('Modality', '?')):3}  {thick:4.1f}mm"
              f"  {ds.get('SeriesDescription', '')}{hint}")
        print(f"  {'':24}  {fs[0].parent.relative_to(root)}")
    print()


def load_series(paths):
    """Read one DICOM series into a volume + a voxel->patient (mm, LPS) affine."""
    slices = [pydicom.dcmread(f) for f in paths]
    if len(slices) < 2:
        sys.exit(f"need at least 2 image slices, found {len(slices)}")
    if len({(ds.Rows, ds.Columns) for ds in slices}) > 1:
        sys.exit("slices have different image sizes; not a single stack")

    iop = np.array(slices[0].ImageOrientationPatient, float)
    if any(np.abs(np.array(ds.ImageOrientationPatient, float) - iop).max() > 1e-3 for ds in slices):
        sys.exit("slices have mixed orientations (a localizer?); pick another series")
    row_dir, col_dir = iop[:3], iop[3:]          # +column index, +row index
    normal = np.cross(row_dir, col_dir)

    # Sort by physical position along the slice normal - not filename, not
    # InstanceNumber, neither of which is guaranteed to be spatial.
    slices.sort(key=lambda ds: np.dot(normal, np.array(ds.ImagePositionPatient, float)))
    pos = np.array([ds.ImagePositionPatient for ds in slices], float)
    gaps = np.diff(pos @ normal)
    if gaps.min() < 1e-3:
        sys.exit("series has overlapping slices (multi-echo / multi-phase?); not supported")
    if np.ptp(gaps) > 0.01 * gaps.mean():
        log(f"warning: uneven slice spacing ({gaps.min():.3f}-{gaps.max():.3f} mm); using mean")
    dk = gaps.mean()
    dr, dc = (float(x) for x in slices[0].PixelSpacing)  # (between rows, between cols)

    vol = np.stack([
        ds.pixel_array.astype(np.float32) * float(getattr(ds, "RescaleSlope", 1))
        + float(getattr(ds, "RescaleIntercept", 0))
        for ds in slices
    ])  # indexed [k, r, c]

    # Columns of A map (k, r, c) voxel indices -> LPS millimetres.
    # Using the step between first and last slice (not the normal) handles
    # gantry-tilted stacks, where slices shear sideways as they advance.
    k_step = (pos[-1] - pos[0]) / (len(slices) - 1)
    A = np.column_stack([k_step, col_dir * dr, row_dir * dc])
    origin = pos[0]

    meta = slices[0]
    log(f"loaded {vol.shape} {getattr(meta, 'Modality', '?')} "
        f"'{getattr(meta, 'SeriesDescription', '')}'  voxel={dk:.3f}x{dr:.3f}x{dc:.3f} mm  "
        f"range=[{vol.min():.0f}, {vol.max():.0f}]")
    return vol, A, origin


# ── 2. segment ────────────────────────────────────────────────────────────────

def body_mask(vol, body_threshold):
    """The patient = largest blob denser than air. Excludes the table/cradle,
    which is separated from the body by foam padding (≈ air on CT)."""
    lab, _ = ndimage.label(ndimage.gaussian_filter(vol, 1) > body_threshold)
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    return lab == sizes.argmax()


def segment(vol, spacing, threshold, sigma_mm, min_island_mm3, max_hole_mm3, body_threshold):
    """Returns a cleaned *continuous* field to contour at `threshold`.

    Contouring the real intensities (not a 0/1 mask) lets marching cubes place
    vertices at sub-voxel positions via partial-volume values, which removes
    the slice-by-slice terracing a binary mask produces.
    """
    voxel_mm3 = float(np.prod(spacing))
    body = body_mask(vol, body_threshold) if body_threshold is not None else None
    if sigma_mm > 0:
        vol = ndimage.gaussian_filter(vol, sigma_mm / np.asarray(spacing))
    mask = vol > threshold
    log(f"threshold > {threshold:g}: {mask.mean() * 100:.2f}% of voxels")

    # Drop floating islands (noise) and anything outside the body (table),
    # but keep every sizeable piece - a skull is several components
    # (cranium, mandible, vertebrae), so "largest only" would be wrong.
    lab, n = ndimage.label(mask)
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    keep = sizes * voxel_mm3 >= min_island_mm3
    if body is not None:
        inside = np.bincount(lab[body].ravel(), minlength=n + 1)
        outside = keep & (inside < 0.5 * sizes)
        keep &= ~outside
        log(f"body mask > {body_threshold:g}: rejected {outside.sum()} piece(s) outside the patient")
    mask = keep[lab]
    log(f"islands: kept {keep.sum()} of {n} (>= {min_island_mm3:g} mm^3)")

    # Fill small enclosed air pockets (trabecular bone, mastoid cells) that
    # would otherwise become thousands of tiny internal shells.
    lab, n = ndimage.label(~mask)
    sizes = np.bincount(lab.ravel())
    small = sizes * voxel_mm3 < max_hole_mm3
    small[0] = False
    small[lab[0, 0, 0]] = False                   # never fill the outside
    holes = small[lab]
    log(f"holes: filled {small.sum()} of {n} (< {max_hole_mm3:g} mm^3)")

    # Keep real intensities everywhere so boundary voxels carry their
    # partial-volume information, but turn rejected pieces into air and
    # filled pockets into solid, so only what we kept can produce a surface.
    lo, hi = float(vol.min()), threshold + abs(threshold) + 1
    field = vol.astype(np.float32, copy=True)
    field[(field > threshold) & ~mask] = lo
    field[holes] = np.maximum(field[holes], hi)
    return field, lo


# ── 3. surface ────────────────────────────────────────────────────────────────

def extract_surface(field, level, background, A, origin):
    # Pad with background so the mesh closes where the scan's field of view
    # cuts through bone - otherwise the top/bottom are open and unprintable.
    padded = np.pad(field, 1, constant_values=background)
    verts, faces, _, _ = measure.marching_cubes(padded, level=level, allow_degenerate=False)
    verts -= 1                                    # undo the pad offset
    mesh = trimesh.Trimesh(verts @ A.T + origin, faces, process=True)
    # Merging coincident vertices can still leave zero-area faces. A vertex
    # listed as its own neighbour corrupts the smoothing Laplacian (rows sum
    # to <1), which drags it toward the scanner origin as a long spike.
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    # Winding depends on skimage's convention *and* the affine's handedness;
    # rather than reason about both, ask the closed mesh which way it faces.
    if mesh.volume < 0:
        mesh.invert()
    log(f"marching cubes: {len(mesh.faces):,} faces")
    return mesh


# ── 4. cleanup ────────────────────────────────────────────────────────────────

def cleanup(mesh, smooth_iters, target_faces):
    if smooth_iters:
        # Taubin alternates shrink/inflate steps, so it removes voxel
        # stair-stepping without the volume loss of plain Laplacian smoothing.
        trimesh.smoothing.filter_taubin(mesh, iterations=smooth_iters)
        log(f"taubin smoothing x{smooth_iters}")
    if target_faces and len(mesh.faces) > target_faces:
        # MeshLab's collapse with preservetopology refuses any edge collapse
        # that would make the mesh non-manifold, so it stays printable.
        ms = pymeshlab.MeshSet()
        ms.add_mesh(pymeshlab.Mesh(mesh.vertices, mesh.faces))
        ms.meshing_decimation_quadric_edge_collapse(
            targetfacenum=target_faces, preservetopology=True, preservenormal=True,
            qualitythr=0.3, optimalplacement=True, planarquadric=True)
        m = ms.current_mesh()
        mesh = trimesh.Trimesh(m.vertex_matrix(), m.face_matrix())
        log(f"decimated -> {len(mesh.faces):,} faces")
    return mesh


# ── 5. export ─────────────────────────────────────────────────────────────────

# LPS (x=Left, y=Posterior, z=Superior) → glTF (Y-up, +Z = front of asset).
LPS_TO_GLTF = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]], float)
FORMATS = (".stl", ".glb", ".obj", ".ply")


def export(mesh, path: Path):
    m = mesh.copy()
    m.apply_translation(-m.bounds.mean(axis=0))   # centre on origin
    if path.suffix.lower() == ".glb":
        T = np.eye(4)
        T[:3, :3] = LPS_TO_GLTF * 0.001           # glTF units are metres
        m.apply_transform(T)
    # STL/OBJ/PLY: leave in LPS millimetres, which is already Z-up.
    path.parent.mkdir(parents=True, exist_ok=True)
    m.export(path)
    log(f"wrote {path}  ({path.stat().st_size / 1e6:.1f} MB)")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=Path,
                   help="folder with DICOM files (a copied disc, or one series' folder)")
    p.add_argument("-s", "--series", type=int,
                   help="which series to convert, by its number in the listing")
    p.add_argument("-o", "--out", type=Path, action="append",
                   help=f"output file, repeatable; format from extension ({' '.join(FORMATS)})")
    p.add_argument("-t", "--threshold", type=float, default=300,
                   help="iso-value; for CT this is Hounsfield units (default 300 = bone)")
    p.add_argument("--sigma", type=float, default=0.4,
                   help="gaussian pre-smoothing in mm, 0 = off (default 0.4)")
    p.add_argument("--min-island", type=float, default=500,
                   help="drop pieces smaller than this, mm^3 (default 500)")
    p.add_argument("--max-hole", type=float, default=200,
                   help="fill cavities smaller than this, mm^3 (default 200)")
    p.add_argument("--smooth", type=int, default=10, help="taubin iterations (default 10)")
    p.add_argument("--faces", type=int, default=600_000,
                   help="decimate to about N faces, 0 = off (default 600000)")
    p.add_argument("--body", type=float, default=-500,
                   help="body-mask threshold to exclude the scanner table (CT HU, default -500)")
    p.add_argument("--no-body", action="store_true", help="disable the body mask (e.g. for MRI)")
    a = p.parse_args()

    if not a.input.is_dir():
        sys.exit(f"not a folder: {a.input}")
    for out in a.out or []:
        if out.suffix.lower() not in FORMATS:
            sys.exit(f"unsupported output format '{out.suffix}'; use one of: {' '.join(FORMATS)}")

    series = find_series(a.input)
    if not series:
        sys.exit(f"no DICOM image series found under {a.input}")
    if not a.out:                                 # no output requested: just list
        print_series(series, a.input)
        print("To convert one:  uv run dicom2mesh.py <folder> -s <#> -o model.stl")
        return
    if a.series is None and len(series) > 1:
        print_series(series, a.input)
        sys.exit("More than one series here: add -s <#> to pick one.")
    if a.series is not None and not 1 <= a.series <= len(series):
        sys.exit(f"--series must be between 1 and {len(series)}")

    _, paths = series[(a.series or 1) - 1]
    vol, A, origin = load_series(paths)
    spacing = np.linalg.norm(A, axis=0)           # mm per voxel along (k, r, c)
    field, bg = segment(vol, spacing, a.threshold, a.sigma, a.min_island, a.max_hole,
                        None if a.no_body else a.body)
    del vol
    mesh = extract_surface(field, a.threshold, bg, A, origin)
    del field
    mesh = cleanup(mesh, a.smooth, a.faces)

    bodies = len(np.unique(trimesh.graph.connected_component_labels(
        mesh.face_adjacency, node_count=len(mesh.faces))))
    log(f"final: {len(mesh.faces):,} faces, {bodies} bodies, watertight={mesh.is_watertight}, "
        f"volume={mesh.volume / 1000:.1f} cm^3, longest edge={mesh.edges_unique_length.max():.1f} mm, "
        f"extent={np.ptp(mesh.vertices, axis=0).round(1)} mm (L-R, A-P, S-I)")
    for out in a.out:
        export(mesh, out)


if __name__ == "__main__":
    main()
