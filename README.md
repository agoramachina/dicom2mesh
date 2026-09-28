# dicom2mesh

Turn the DICOM files from a CT (or MRI) disc into a 3D model you can view, share, or 3D-print.

```
DICOM slices ──► 3D volume ──► segmentation ──► surface mesh ──► cleanup ──► STL / GLB / OBJ / PLY
```

One self-contained Python script. It works best for **bone from CT** (e.g. a skull from a head CT). MRI works only roughly for now; see [MRI](#mri) below.

> **Not a medical device.** This is a hobby / proof-of-concept tool. Don't use it for diagnosis, surgical planning, or anything clinical.

---

## Install

You only need [**uv**](https://docs.astral.sh/uv/). It fetches Python and every dependency automatically the first time you run the script (a few hundred MB, cached afterwards). No virtualenv, no `pip install`, and no compiler.

| OS | Install uv |
|---|---|
| **Windows** (PowerShell) | `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 \| iex"` or `winget install --id=astral-sh.uv -e` |
| **macOS** | `curl -LsSf https://astral.sh/uv/install.sh \| sh` or `brew install uv` |
| **Linux** | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |

Then download `dicom2mesh.py` (and optionally `preview.py`) into a folder and open a terminal there.

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

MRI brightness has **no fixed scale**. It changes with the scanner, the sequence, and even position within the scan, and bone is *dark*. So a single threshold only gives a rough result. If you want to experiment, use `--no-body` and choose a threshold from the `range=[...]` printed when the series loads. A proper **brain** model needs real tissue segmentation, which this script doesn't do (yet).

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

- Threshold-based: great for CT bone, rough for MRI; no tissue classification.
- Very thin bone (e.g. eye-socket walls) can come out with holes, which is normal for CT skull models.
- Enhanced multi-frame DICOM isn't supported yet.
- Tested end-to-end on Linux with one head CT (Siemens) and its reformats. Dependencies are verified to install from prebuilt wheels on Windows x64, macOS (Intel & Apple Silicon), and Linux x64/ARM for Python 3.10–3.14, but it hasn't yet been run on Windows or macOS hardware. Reports welcome.

## How this was made

This was written with [Claude Code](https://claude.com/claude-code) (Opus 5.5) as an experiment in what an AI assistant can and can't do for this kind of task, starting from a real set of scans from a patient CD.
