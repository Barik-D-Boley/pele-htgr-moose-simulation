"""
Render temperature and von Mises stress cross-sections at the hottest height
for both MOOSE cases, side by side, with ParaView.

Run with ParaView's Python (not plain python), from the repository root:
    pvpython post-processing/render_fields.py
    (Windows: "C:/Program Files/ParaView 6.1.1/bin/pvpython.exe" post-processing/render_fields.py)

Writes images/field_temperature.png and images/field_vonmises.png.
Both cases share one color range per field, so colors compare directly.
"""
import os
import sys

from paraview.simple import (AssignViewToLayout, Calculator, ColorBy, CreateLayout, CreateRenderView,
                             ExodusIIReader, ExtractBlock, FeatureEdges, GetColorTransferFunction,
                             GetScalarBar, MergeBlocks, Render, SaveScreenshot, Show, Slice, Text)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = [('uniform', 'A: uniform heat'), ('openmc', 'B: OpenMC heat')]
# Hottest fuel is at the bottom of the column (helium outlet). Slice 5 mm up
# so the cut sits inside the first layer of elements, not on the boundary face.
Z_SLICE = 0.005  # m
FIELDS = [
    # display array, point/cell, title, units, image
    ('T_C', 'POINTS', 'Temperature', '°C', 'field_temperature.png'),
    ('vm_MPa', 'CELLS', 'von Mises stress', 'MPa', 'field_vonmises.png'),
]


def load(case):
    path = os.path.join(ROOT, 'moose', f'pele_{case}_out.e')
    if not os.path.exists(path):
        sys.exit(f'Missing {path}. Run the MOOSE cases first.')
    r = ExodusIIReader(FileName=[path])
    r.PointVariables = ['temperature']
    r.ElementVariables = ['vonmises_stress']
    r.UpdatePipeline(time=r.TimestepValues[-1] if r.TimestepValues else 0.0)
    m = MergeBlocks(Input=r)
    s = Slice(Input=m)
    s.SliceType = 'Plane'
    s.SliceType.Origin = [0.0, 0.0, Z_SLICE]
    s.SliceType.Normal = [0.0, 0.0, 1.0]
    # Convert to display units: temperature in C, stress in MPa
    # (one branch per field: chaining a point and a cell Calculator drops the cell array)
    t_end = r.TimestepValues[-1] if r.TimestepValues else 0.0
    t = Calculator(Input=s, AttributeType='Point Data', ResultArrayName='T_C',
                   Function='temperature - 273.15')
    vm = Calculator(Input=s, AttributeType='Cell Data', ResultArrayName='vm_MPa',
                    Function='vonmises_stress * 1e-6')
    t.UpdatePipeline(time=t_end)
    vm.UpdatePipeline(time=t_end)
    # Outline of the fuel compact, drawn over the field
    fb = ExtractBlock(Input=r, Selectors=['/Root/ElementBlocks/fuel_compact'])
    fs = Slice(Input=MergeBlocks(Input=fb))
    fs.SliceType = 'Plane'
    fs.SliceType.Origin = [0.0, 0.0, Z_SLICE]
    fs.SliceType.Normal = [0.0, 0.0, 1.0]
    edges = FeatureEdges(Input=fs)
    edges.UpdatePipeline(time=t_end)
    return {'T_C': t, 'vm_MPa': vm, 't_end': t_end, 'fuel_edge': edges}


def data_range(src, name, assoc):
    info = src.GetPointDataInformation() if assoc == 'POINTS' else src.GetCellDataInformation()
    return info.GetArray(name).GetRange()


def main():
    slices = {case: load(case) for case, _ in CASES}
    out_dir = os.path.join(ROOT, 'images')
    for name, assoc, title, unit, png in FIELDS:
        lo = min(data_range(src[name], name, assoc)[0] for src in slices.values())
        hi = max(data_range(src[name], name, assoc)[1] for src in slices.values())
        layout = CreateLayout(title)
        views = []
        for i, (case, label) in enumerate(CASES):
            v = CreateRenderView()
            v.ViewSize = [700, 420]
            v.ViewTime = slices[case]['t_end']   # final (solved) step, not the initial guess
            v.OrientationAxesVisibility = 0
            v.Background = [1, 1, 1]
            v.UseColorPaletteForBackground = 0
            d = Show(slices[case][name], v)
            ColorBy(d, (assoc, name))
            lut = GetColorTransferFunction(name)
            # Preset names differ between ParaView versions; use the first one found
            presets = (['Inferno (matplotlib)', 'Inferno'] if name == 'T_C'
                       else ['Viridis (matplotlib)', 'Viridis', 'viridis'])
            for preset in presets:
                try:
                    lut.ApplyPreset(preset, True)
                    break
                except RuntimeError:
                    continue
            lut.RescaleTransferFunction(lo, hi)
            bar = GetScalarBar(lut, v)
            bar.Title = f'{title} ({unit})'
            bar.ComponentTitle = ''
            bar.TitleColor = bar.LabelColor = [0.2, 0.2, 0.2]
            bar.RangeLabelFormat = '%-.0f'
            bar.LabelFormat = '%-.0f'
            bar.WindowLocation = 'Any Location'
            bar.Position = [0.88, 0.12]
            bar.ScalarBarLength = 0.75
            bar.Visibility = 1
            e = Show(slices[case]['fuel_edge'], v)
            e.ColorArrayName = ['POINTS', '']
            e.AmbientColor = e.DiffuseColor = [1, 1, 1]
            e.LineWidth = 2.0
            t = Text(Text=f'{label}   (z = {Z_SLICE * 100:.1f} cm, hottest height; white line = fuel compact)')
            td = Show(t, v)
            td.Color = [0.2, 0.2, 0.2]
            td.FontSize = 14
            td.WindowLocation = 'Upper Left Corner'
            views.append(v)
        layout.SplitHorizontal(0, 0.5)
        AssignViewToLayout(view=views[0], layout=layout, hint=1)
        AssignViewToLayout(view=views[1], layout=layout, hint=2)
        for v in views:
            Render(v)   # first render resets the camera, so set ours after it
            # Same camera in both panels; leave room on the right for the color bar
            v.CameraPosition = [0.0123, 0.0047, 0.2]
            v.CameraFocalPoint = [0.0123, 0.0047, Z_SLICE]
            v.CameraViewUp = [0, 1, 0]
            v.CameraParallelProjection = 1
            v.CameraParallelScale = 0.0076
            Render(v)
        path = os.path.join(out_dir, png)
        SaveScreenshot(path, layout, ImageResolution=[1400, 420])
        print(f'Saved images/{png}  ({title}: {lo:.1f} to {hi:.1f} {unit})')


if __name__ == '__main__':
    main()
