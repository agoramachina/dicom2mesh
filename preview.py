# /// script
# dependencies = ["vtk", "numpy", "pillow"]
# ///
"""Offscreen turntable-style preview of an STL (assumes LPS mm, Z-up)."""
import sys, numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
from PIL import Image
src, out = sys.argv[1], sys.argv[2]
r = vtk.vtkSTLReader(); r.SetFileName(src); r.Update()
c = np.array(r.GetOutput().GetCenter())
views = {"front": (0,-1,0), "right": (-1,0,0), "3/4": (-0.8,-1,0.35), "top": (0,-0.15,1), "back": (0.4,1,0.2), "below": (0.2,-0.4,-1)}
tiles = []
for name,(dx,dy,dz) in views.items():
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(r.GetOutputPort())
    a = vtk.vtkActor(); a.SetMapper(m)
    p = a.GetProperty(); p.SetColor(0.93,0.89,0.80); p.SetSpecular(0.25); p.SetSpecularPower(20); p.SetAmbient(0.12)
    ren = vtk.vtkRenderer(); ren.AddActor(a); ren.SetBackground(0.08,0.09,0.11); ren.SetBackground2(0.22,0.24,0.28); ren.GradientBackgroundOn()
    win = vtk.vtkRenderWindow(); win.SetOffScreenRendering(1); win.AddRenderer(ren); win.SetSize(600,600); win.SetMultiSamples(8)
    cam = ren.GetActiveCamera(); d = np.array([dx,dy,dz],float); d/=np.linalg.norm(d)
    cam.SetFocalPoint(*c); cam.SetPosition(*(c+d*500)); cam.SetViewUp(0,1,0) if abs(dz)>0.9 else cam.SetViewUp(0,0,1)
    ren.ResetCamera(); cam.Zoom(1.25); ren.ResetCameraClippingRange()
    win.Render()
    w = vtk.vtkWindowToImageFilter(); w.SetInput(win); w.Update()
    img = w.GetOutput(); h, wd, _ = img.GetDimensions()
    arr = vtk_to_numpy(img.GetPointData().GetScalars()).reshape(wd, h, -1)[::-1, :, :3]
    tiles.append(Image.fromarray(np.ascontiguousarray(arr)))
W = Image.new("RGB", (1800, 1200))
for i,t in enumerate(tiles): W.paste(t, ((i%3)*600, (i//3)*600))
W.save(out); print("saved", out)
