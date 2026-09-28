# dicom2mesh

Turn the DICOM files from a CT (or MRI) disc into a 3D model you can view, share, or 3D-print.

```
DICOM slices ──► 3D volume ──► segmentation ──► surface mesh ──► cleanup ──► STL / GLB / OBJ / PLY
```

Two scripts:

- **`dicom2mesh.py`**: **bone from CT** (e.g. a skull from a head CT). Threshold-based, fast, no ML.
- **`brain2mesh.py`**: **a brain surface from a T1 MRI** (e.g. MPRAGE). Uses small ML models to find the brain; see [Brain from MRI](#brain-from-mri).

> **Not a medical device.** This is a hobby / proof-of-concept tool. Don't use it for diagnosis, surgical planning, or anything clinical.

---

## Install

You only need [**uv**](https://docs.astral.sh/uv/). It fetches Python and every dependency automatically the first time you run the script (a few hundred MB, cached afterwards). No virtualenv, no `pip install`, and no compiler.

| OS | Install uv |
|---|---|
| **Windows** (PowerShell) | `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 \| iex"` or `winget install --id=astral-sh.uv -e` |
| **macOS** | `curl -LsSf https://astral.sh/uv/install.sh \| sh` or `brew install uv` |
| **Linux** | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |

Then download the scripts into **one folder** (`brain2mesh.py` reuses `dicom2mesh.py`, so keep them together) and open a terminal there.

---

## Quick start

### 1. Copy your disc to a folder

Copy **the whole disc** to your computer (e.g. `C:\scans\head-ct` or `~/scans/head-ct`). Reading straight off a CD is slow, and a bad copy is the most common problem. The script warns you if it finds empty files, which usually means the copy didn't finish.

### 2. List what's on it

```
uv run dicom2mesh.py C:\scans\head-ct
```

A typical disc holds several *series* (separate scans or reconstructions):

```
Found 5 image series under C:\scans\head-ct:

   #  slices  mod   slice  description / folder
   1       1  CT    0.6mm  Topogram  0.6  T20f
                            DICOM/0001/0002
   2      32  CT    5.0mm  AXIAL_WO_STD
                            DICOM/0001/0003
   3     157  CT    1.0mm  AXIAL_WO_BONE   <- good for 3D
                            DICOM/0001/0004
   4     100  CT    2.0mm  COR
                            DICOM/0001/0005
   5      85  CT    2.0mm  SAG
                            DICOM/0001/0006
```

### 3. Convert one

```
uv run dicom2mesh.py C:\scans\head-ct -s 3 -o skull.stl -o skull.glb
```

This takes about a minute and prints a line per step. The last line tells you whether the result is good:

```
final: 600,000 faces, 3 bodies, watertight=True, volume=412.7 cm^3, longest edge=6.3 mm, extent=[136.5 178. 145.5] mm (L-R, A-P, S-I)
```

`watertight=True` means it's ready to print. The `extent` is the real-world size in millimetres; an adult skull is about 140 × 180 mm.

### 4. (Optional) Preview it

```
uv run preview.py skull.stl
```

This writes `skull_preview.png`, with front, side, three-quarter, top, back, and underside views. It renders offscreen, so no window opens.

---

## Choosing a series

- **Thin slices matter most.** Pick the series with the most slices and the thinnest slice thickness (`<- good for 3D` marks series with 50+ slices at ≤ 1.5 mm). Thick slices show up as visible terraces in the model.
- **CT reconstruction kernel.** Head CTs often come in a *bone* version (sharp: names like `BONE`, `H60`, `B70`) and a *soft-tissue* version (smooth: `STD`, `J40`, `H30`). Use the **bone** one for skulls.
- **COR / SAG** series are usually reformats of the same scan: they work, but typically have thicker slices than the original axial series.
- **Localizers / scouts** (1–3 slices) can't be converted.

---

## Output formats

The format comes from the file extension, and you can pass `-o` several times.

| Format | Best for | Units / orientation |
|---|---|---|
| **`.stl`** | 3D printing. Every slicer reads it. | millimetres, Z-up |
| **`.glb`** | Viewing & sharing: Blender, Windows 3D Viewer, web viewers (e.g. drag onto [gltf-viewer.donmccurdy.com](https://gltf-viewer.donmccurdy.com/)) | metres, Y-up (glTF standard) |
| **`.obj`** | Importing into other 3D software | millimetres, Z-up |
| **`.ply`** | Compact, for mesh tools (MeshLab, CloudCompare) | millimetres, Z-up |

All outputs are centred on the origin, at real-world (1:1) scale.

---

## Options

| Option | Default | What it does |
|---|---|---|
| `-s N`, `--series N` | | Which series to convert (number from the listing) |
| `-o FILE` | | Output file; repeatable |
| `-t`, `--threshold` | `300` | Surface level. For CT this is in Hounsfield units (HU) |
| `--sigma` | `0.4` | Pre-smoothing in mm. Higher = smoother, less detail. `0` = off |
| `--smooth` | `10` | Mesh smoothing iterations (Taubin; doesn't shrink the model) |
| `--faces` | `600000` | Reduce to about this many triangles. `0` = keep all (millions) |
| `--min-island` | `500` | Discard loose pieces smaller than this (mm³) |
| `--max-hole` | `200` | Fill enclosed cavities smaller than this (mm³) |
| `--body` | `-500` | Threshold for the "body" mask that removes the scanner table (CT) |
| `--no-body` | | Turn the body mask off (use for MRI) |

### CT thresholds (starting points)

CT values are calibrated: air is −1000 HU and water is 0 HU on every scanner. That's why a single threshold works.

| Tissue | Approx. HU | Try |
|---|---|---|
| Bone | 300 to 2000+ | `-t 300` (default) |
| Dense bone only | 700+ | `-t 700` (thinner, cleaner, but more holes) |
| Skin surface | soft tissue ≈ 0–80, fat ≈ −100 | `-t -300` ⚠️ produces a recognizable face; see [Privacy](#privacy) |

### MRI

MRI brightness has **no fixed scale**. It changes with the scanner, the sequence, and even position within the scan, and bone is *dark*. So a single threshold only gives a rough result. If you want to experiment, use `--no-body` and choose a threshold from the `range=[...]` printed when the series loads. **For a brain, use [`brain2mesh.py`](#brain-from-mri) instead.**

---

## Brain from MRI

`brain2mesh.py` makes a **brain surface with its folds** (the *pial* surface: cerebrum, cerebellum, and brainstem) from a **3D T1-weighted MRI**. On a disc, that's usually a series with `MPRAGE`, `T1 3D`, `SPGR`, or `BRAVO` in its name, about 1 mm thick. Scans with contrast (gadolinium) work too.

```
uv run brain2mesh.py C:\scans\head-mri                       # list the series
uv run brain2mesh.py C:\scans\head-mri -s 10 -o brain.stl -o brain.glb
```

How it works:

1. **Skull stripping** and **tissue labelling** with two [brainchop](https://github.com/neuroneural/brainchop-cli) models (`mindgrab`, then `subcortical`). These are small MeshNet networks, trained on synthetic scans so they don't depend on one particular contrast. They download automatically from GitHub on first use.
2. **Placing the surface:** the labels decide *what* is brain; the T1 brightness decides *exactly where* its edge is. The script contours halfway between grey-matter and CSF brightness, measured on your own scan, so the folds come out crisp instead of voxel-lumpy. The surface is only allowed to move 1 voxel outward from the labels, so blood vessels and the brain's outer membrane don't get stuck to it.
3. The same cleanup, smoothing, and export steps as `dicom2mesh.py`.

| Option | Default | What it does |
|---|---|---|
| `--sigma` | `0.7` | Pre-smoothing of the T1 in mm. Higher = smoother folds |
| `--smooth` | `15` | Mesh smoothing iterations |
| `--reach` | `1` | Voxels the surface may move outward from the labels (higher = fuller, but more debris) |
| `--device` | auto | Force brainchop's compute device, e.g. `--device CPU` |
| `--keep DIR` | | Keep the intermediate NIfTI files (stripped scan, labels) |

**Hardware: you effectively need a GPU.** brainchop picks one automatically. On an NVIDIA laptop GPU with ~3 GB free, both models took about 30 seconds in total. If you see `out of GPU memory`, close other GPU-heavy apps (3D viewers, slicers, games) and try again. `--device CPU` exists, but in testing (brainchop 0.2.5) the CPU run hadn't finished even the *first* model after 20 minutes, running on a single core.

---

## 3D printing tips

- Open the `.stl` in your slicer (PrusaSlicer, Bambu Studio, OrcaSlicer, Cura, …). It's already in millimetres at life size, so scale it down to save time and filament.
- A skull has overhangs everywhere (eye sockets, skull base), so turn on supports; *tree/organic* supports come off most cleanly.
- The model is watertight (closed), which is what slicers need. The flat surface at the bottom is where the scan ended.

---

## Privacy

- **DICOM files contain personal data**: name, date of birth, patient ID, the facility, dates. **Never upload or commit your scan files.**
- The meshes this script writes contain **no metadata**, only geometry.
- But geometry can identify you as well: a **skin-surface** model (`-t -300`) is your face. Skull models are much less identifying, but still yours to decide about.
- Keep scans and generated models **outside** this folder. The included `.gitignore` also blocks common mesh and DICOM file types as a backstop, but many discs use extensionless file names that no pattern can catch.

---

## Troubleshooting

| Message / problem | Likely cause & fix |
|---|---|
| `skipped N empty (0-byte) files` | The copy from the disc is incomplete. Copy it again (another drive, or a recovery tool like `ddrescue` for scratched discs). |
| `no DICOM image series found` | The disc may use *enhanced multi-frame* DICOM (one big file per series), which isn't supported yet, or the folder isn't the disc copy. |
| `slices have mixed orientations` | That series is a localizer; pick another. |
| `overlapping slices` | A multi-echo / multi-phase series; pick another. |
| `uneven slice spacing` (warning) | Usually harmless; if the model looks sheared or stretched, pick another series. |
| Runs out of memory | A 157-slice 512×512 CT peaks at about 2 GB of RAM; usage scales with slice count. Close other apps. |
| **Windows:** path errors | Put paths with spaces in quotes, but **don't end a quoted path with a backslash**: `"C:\my scans\"` breaks, because `\"` is read as an escaped quote. Use `"C:\my scans"` or forward slashes. |
| **Linux:** dependency won't install | Needs a reasonably current distro (glibc 2.31+, e.g. Ubuntu 20.04+). |

Compressed discs (JPEG Lossless, JPEG 2000, JPEG-LS, RLE) are supported via GDCM.

---

## How it works

1. **Find**: recursively reads DICOM *headers only* and groups files by `SeriesInstanceUID`.
2. **Load**: sorts slices by physical position (`ImagePositionPatient`) rather than by filename, converts pixel values to real units (HU for CT), and builds a voxel→millimetre transform. This handles oblique and gantry-tilted scans.
3. **Segment**: thresholds, keeps only pieces inside the patient's body (dropping the scanner table), removes small specks, and fills small internal air pockets.
4. **Surface**: marching cubes runs on the *actual intensity values* rather than a yes/no mask, so vertices land at sub-voxel positions. This greatly reduces stair-stepping. The volume is padded so the mesh closes where the scan ends.
5. **Clean up**: removes degenerate triangles, applies Taubin smoothing, reduces triangle count with MeshLab's topology-preserving decimation (keeps the mesh printable), and fixes inside-out normals using the sign of the enclosed volume.
6. **Export**: centres the model and converts DICOM patient coordinates (LPS) to each format's conventions.

## Limitations

- `dicom2mesh.py` is threshold-based: great for CT bone, rough for MRI.
- Very thin bone (e.g. eye-socket walls) can come out with holes, which is normal for CT skull models.
- `brain2mesh.py` makes one surface for the whole brain: the hemispheres aren't split, and there's no separate white-matter surface yet. At 1 mm resolution, the tightest folds can come out bridged. It's a faithful sculpture of your brain, not a research-grade reconstruction like FreeSurfer.
- brainchop's `subcortical` model ships without label names; the two this script relies on (cortex, ventricles) were identified by inspecting their size, position, and brightness.
- Enhanced multi-frame DICOM isn't supported yet.
- Tested end-to-end on Linux with one head CT (Siemens) and its reformats, and one contrast-enhanced 3 T MPRAGE (Siemens) on an NVIDIA GPU. Dependencies are verified to install from prebuilt wheels on Windows x64, macOS (Intel & Apple Silicon), and Linux x64/ARM for Python 3.10–3.14, but neither script has been run on Windows or macOS hardware yet. Reports welcome.

## Credits

Brain segmentation uses [brainchop](https://github.com/neuroneural/brainchop-cli) (MIT) and its MeshNet models, from the neuroneural group.

## How this was made

This was written with [Claude Code](https://claude.com/claude-code) (Opus 5.5) as an experiment in what an AI assistant can and can't do for this kind of task, starting from a real set of scans from a patient CD.
