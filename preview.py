# /// script
# requires-python = ">=3.10"
# dependencies = ["vtk", "numpy", "pillow"]
# ///
"""
preview - render six views of a dicom2mesh output into one PNG.

    uv run preview.py skull.stl                 # writes skull_preview.png
    uv run preview.py skull.stl views.png

Takes the .stl/.ply/.obj that dicom2mesh writes (patient coordinates in mm,
Z = up). Renders offscreen, so no window opens.
"""

import sys
from pathlib import Path

import numpy as np
import vtk
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy

READERS = {".stl": vtk.vtkSTLReader, ".ply": vtk.vtkPLYReader, ".obj": vtk.vtkOBJReader}

# Camera directions in patient coordinates (x=Left, y=Posterior, z=Superior).
VIEWS = {
    "front": (0, -1, 0), "right": (-1, 0, 0), "front-right": (-0.8, -1, 0.35),
    "top": (0, -0.15, 1), "back": (0.4, 1, 0.2), "below": (0.2, -0.4, -1),
}
TILE = 600


def render(src: Path, out: Path):
    reader = READERS[src.suffix.lower()]()
    reader.SetFileName(str(src))
    reader.Update()
    centre = np.array(reader.GetOutput().GetCenter())

    tiles = []
    for name, d in VIEWS.items():
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(reader.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        prop = actor.GetProperty()
        prop.SetColor(0.93, 0.89, 0.80)          # bone
        prop.SetSpecular(0.25)
        prop.SetSpecularPower(20)
        prop.SetAmbient(0.12)

        label = vtk.vtkTextActor()
        label.SetInput(name)
        label.SetPosition(12, 12)
        label.GetTextProperty().SetFontSize(22)
        label.GetTextProperty().SetColor(0.75, 0.78, 0.82)

        ren = vtk.vtkRenderer()
        ren.AddActor(actor)
        ren.AddViewProp(label)
        ren.SetBackground(0.08, 0.09, 0.11)
        ren.SetBackground2(0.22, 0.24, 0.28)
        ren.GradientBackgroundOn()

        win = vtk.vtkRenderWindow()
        win.SetOffScreenRendering(1)
        win.AddRenderer(ren)
        win.SetSize(TILE, TILE)
        win.SetMultiSamples(8)

        d = np.array(d, float) / np.linalg.norm(d)
        cam = ren.GetActiveCamera()
        cam.SetFocalPoint(*centre)
        cam.SetPosition(*(centre + d * 500))
        cam.SetViewUp(0, 1, 0) if abs(d[2]) > 0.9 else cam.SetViewUp(0, 0, 1)
        ren.ResetCamera()
        cam.Zoom(1.25)
        ren.ResetCameraClippingRange()
        win.Render()

        grab = vtk.vtkWindowToImageFilter()
        grab.SetInput(win)
        grab.Update()
        img = grab.GetOutput()
        w, h, _ = img.GetDimensions()
        rgb = vtk_to_numpy(img.GetPointData().GetScalars()).reshape(h, w, -1)[::-1, :, :3]
        tiles.append(Image.fromarray(np.ascontiguousarray(rgb)))
        win.Finalize()

    sheet = Image.new("RGB", (3 * TILE, 2 * TILE))
    for i, t in enumerate(tiles):
        sheet.paste(t, ((i % 3) * TILE, (i // 3) * TILE))
    sheet.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3) or Path(sys.argv[1]).suffix.lower() not in READERS:
        sys.exit(__doc__)
    src = Path(sys.argv[1])
    render(src, Path(sys.argv[2]) if len(sys.argv) == 3 else src.with_name(src.stem + "_preview.png"))
