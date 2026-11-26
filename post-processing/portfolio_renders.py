"""
Presentation renders of the MOOSE results (portfolio / slides), 1920 x 1080.

The model is a quarter of a unit cell. For display it is mirrored about
x = 0 and y = 0 to show the full cell: one coolant hole in the middle and a
quarter fuel compact in each corner.

Run with ParaView's Python, from the repository root:
    pvpython post-processing/portfolio_renders.py
    (Windows: "C:/Program Files/ParaView 6.1.1/bin/pvpython.exe" post-processing/portfolio_renders.py)

Writes images/portfolio/*.png.
"""
import os
import sys

from paraview.simple import (Calculator, Clip, ColorBy, CreateRenderView, ExodusIIReader, ExtractBlock,
                             CleantoGrid, ExtractSurface, FeatureEdges,
                             GetColorTransferFunction, GetScalarBar, MergeBlocks, Reflect, Render,
                             SaveScreenshot, Show, Transform, WarpByVector, Delete, Hide)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'images', 'portfolio')
SIZE = [1920, 1080]
SLAB = 0.013         # m of column shown in the close-ups (two element layers at the hottest end)
Z_SCALE = 0.12       # height scaling for the full-column renders (80 cm column -> 9.6 cm)


def load(case):
    path = os.path.join(ROOT, 'moose', f'pele_{case}_out.e')
    if not os.path.exists(path):
        sys.exit(f'Missing {path}. Run the MOOSE cases first.')
    r = ExodusIIReader(FileName=[path])
    r.PointVariables = ['temperature', 'disp_']   # the reader groups disp_x/y/z into one vector
    r.ElementVariables = ['vonmises_stress', 'hydrostatic_stress', 'power_q']
    r.ApplyDisplacements = 0
    t_end = r.TimestepValues[-1]
    r.UpdatePipeline(t_end)
    return r, t_end


def full_cell(src):
    """Mirror the quarter cell about x = 0 and then y = 0."""
    # Mirror about the planes x = 0 and y = 0 exactly ('X Min' would use each
    # dataset's own bounds, which is wrong for the fuel block on its own).
    m = MergeBlocks(Input=src)
    rx = Reflect(Input=m)
    rx.Plane = 'X'
    rx.Center = 0.0
    rx.CopyInput = 1
    ry = Reflect(Input=MergeBlocks(Input=rx))
    ry.Plane = 'Y'
    ry.Center = 0.0
    ry.CopyInput = 1
    return MergeBlocks(Input=ry)


def units(src, fuel_only=False):
    """Display units: temperature in C, stresses in MPa, power in MW/m3, displacement in um.
    power_q exists only in the fuel block, so it is converted only for fuel-only inputs."""
    c = Calculator(Input=src, AttributeType='Point Data', ResultArrayName='T_C', Function='temperature - 273.15')
    c = Calculator(Input=c, AttributeType='Point Data', ResultArrayName='disp_um', Function='mag(disp_) * 1e6')
    cells = [('vm_MPa', 'vonmises_stress * 1e-6'), ('hydro_MPa', 'hydrostatic_stress * 1e-6')]
    if fuel_only:
        cells.append(('q_MW', 'power_q * 1e-6'))
    for name, expr in cells:
        c = Calculator(Input=c, AttributeType='Cell Data', ResultArrayName=name, Function=expr)
    return c


def bottom_slab(src):
    clip = Clip(Input=src)
    clip.ClipType = 'Plane'
    clip.ClipType.Origin = [0, 0, SLAB]
    clip.ClipType.Normal = [0, 0, 1]
    clip.Invert = 1          # keep z < SLAB
    clip.Crinkleclip = 1     # keep whole elements so the mesh edges stay clean
    return clip


def view():
    v = CreateRenderView()
    v.ViewSize = SIZE
    v.OrientationAxesVisibility = 0
    v.UseColorPaletteForBackground = 0
    v.Background = [1, 1, 1]
    return v


def color(d, v, assoc, name, title, preset, fmt='%-.0f', rng=None):
    ColorBy(d, (assoc, name))
    lut = GetColorTransferFunction(name)
    for p in ([preset] if isinstance(preset, str) else preset):
        try:
            lut.ApplyPreset(p, True)
            break
        except RuntimeError:
            continue
    if rng:
        lut.RescaleTransferFunction(*rng)
    else:
        d.RescaleTransferFunctionToDataRange(False, True)
    bar = GetScalarBar(lut, v)
    bar.Title = title
    bar.ComponentTitle = ''
    bar.TitleColor = bar.LabelColor = [0, 0, 0]
    bar.TitleFontSize = 30
    bar.TitleBold = 1
    bar.LabelFontSize = 26
    bar.RangeLabelFormat = fmt
    bar.LabelFormat = fmt
    bar.WindowLocation = 'Any Location'
    bar.Position = [0.80, 0.30]
    bar.ScalarBarLength = 0.40
    bar.ScalarBarThickness = 22
    bar.Visibility = 1
    return lut


def iso_camera(v, center, dist, elev=(1.0, -1.6, 1.15)):
    Render(v)   # the first render resets the camera, so set ours after it
    v.CameraFocalPoint = center
    v.CameraPosition = [center[0] + elev[0] * dist, center[1] + elev[1] * dist, center[2] + elev[2] * dist]
    v.CameraViewUp = [0, 0, 1]
    v.CameraViewAngle = 30


def save(v, name):
    Render(v)
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    SaveScreenshot(path, v, ImageResolution=SIZE)
    print(f'Saved images/portfolio/{name}')


def outline(src, v, color=(0.05, 0.05, 0.05), width=2.0):
    """Draw the sharp edges of a block (e.g. the fuel compacts) as lines."""
    e = FeatureEdges(Input=ExtractSurface(Input=src))
    e.FeatureAngle = 60
    d = Show(e, v)
    d.ColorArrayName = ['POINTS', '']
    d.AmbientColor = d.DiffuseColor = list(color)
    d.LineWidth = width
    return d


def block(r, name):
    return ExtractBlock(Input=r, Selectors=[f'/Root/ElementBlocks/{name}'])


# Camera looking down on the wafer at about 55 degrees, from the front-right
WAFER = dict(center=[0.0055, -0.002, SLAB * 1.6], dist=0.047, elev=(0.55, -0.95, 1.45))


def main():
    r, t_end = load('openmc')
    cell = units(full_cell(r))
    slab = bottom_slab(cell)
    fuel_slab = bottom_slab(full_cell(block(r, 'fuel_compact')))

    # 1-3: the hottest end of the column, full cell, fuel compacts outlined
    shots = [
        ('POINTS', 'T_C', 'Temperature (°C)', ['Cool to Warm (Extended)'], 'pele_paraview_temp.png'),
        ('CELLS', 'vm_MPa', 'von Mises stress (MPa)', ['Turbo', 'Jet', 'Rainbow Uniform'], 'pele_paraview_vonmises.png'),
        ('CELLS', 'hydro_MPa', 'Hydrostatic stress (MPa)', ['Cool to Warm'], 'pele_paraview_hydro.png'),
    ]
    for assoc, name, title, preset, png in shots:
        v = view()
        v.ViewTime = t_end
        d = Show(slab, v)
        lut = color(d, v, assoc, name, title, preset)
        if name == 'hydro_MPa':
            # diverging map centered on zero: blue = compression, red = tension
            rng = slab.GetCellDataInformation().GetArray(name).GetRange()
            m = max(abs(rng[0]), abs(rng[1]))
            lut.RescaleTransferFunction(-m, m)
            GetScalarBar(lut, v).Title = 'Hydrostatic stress (MPa)\n(− compression, + tension)'
        outline(fuel_slab, v)
        iso_camera(v, WAFER['center'], WAFER['dist'], WAFER['elev'])
        save(v, png)
        Delete(v)

    # 4: displacement, exaggerated so the thermal expansion shape is visible
    warp = WarpByVector(Input=slab)
    warp.Vectors = ['POINTS', 'disp_']
    warp.ScaleFactor = 20.0
    v = view()
    v.ViewTime = t_end
    d = Show(warp, v)
    color(d, v, 'POINTS', 'disp_um', 'Displacement (µm)\nshape exaggerated ×20', ['Viridis (matplotlib)', 'Viridis'])
    ghost = Show(slab, v)          # undeformed shape for reference
    ghost.SetRepresentationType('Feature Edges')
    ghost.ColorArrayName = ['POINTS', '']
    ghost.AmbientColor = ghost.DiffuseColor = [0.45, 0.45, 0.45]
    iso_camera(v, WAFER['center'], WAFER['dist'] * 1.15, WAFER['elev'])
    save(v, 'pele_paraview_disp.png')
    Delete(v)

    # 5-6: full 80 cm column, height scaled so the axial profile fits on screen
    tall = Transform(Input=cell)
    tall.Transform.Scale = [1.0, 1.0, Z_SCALE]
    col_center = [0.0, 0.0, 0.8 * Z_SCALE / 2]
    col_cam = dict(dist=0.14, elev=(0.75, -1.25, 0.55))

    v = view()
    v.ViewTime = t_end
    d = Show(tall, v)
    color(d, v, 'POINTS', 'T_C', 'Temperature (°C)\nheight scaled ×0.12', ['Cool to Warm (Extended)'])
    iso_camera(v, col_center, col_cam['dist'], col_cam['elev'])
    save(v, 'pele_paraview_column_temp.png')
    Delete(v)

    v = view()
    v.ViewTime = t_end
    fuel = Transform(Input=units(full_cell(block(r, 'fuel_compact')), fuel_only=True))
    fuel.Transform.Scale = [1.0, 1.0, Z_SCALE]
    graphite = Transform(Input=full_cell(block(r, 'matrix_graphite')))
    graphite.Transform.Scale = [1.0, 1.0, Z_SCALE]
    g = Show(graphite, v)                 # graphite as gray outline edges for context
    g.SetRepresentationType('Feature Edges')
    g.ColorArrayName = ['POINTS', '']
    g.AmbientColor = g.DiffuseColor = [0.55, 0.55, 0.55]
    g.LineWidth = 1.5
    d = Show(fuel, v)
    color(d, v, 'CELLS', 'q_MW', 'Fission power (MW/m³)\nheight scaled ×0.12', ['Inferno (matplotlib)', 'Inferno'])
    iso_camera(v, col_center, col_cam['dist'], col_cam['elev'])
    save(v, 'pele_paraview_power.png')
    Delete(v)

    # Geometry, CAD style: smooth shaded blocks with edge lines, bottom 6 cm
    # of the column (taller than the result close-ups, so it reads as a column)
    v = view()
    v.ViewTime = t_end
    geo_clip = 0.06
    for name, rgb in (('matrix_graphite', [0.62, 0.64, 0.68]), ('fuel_compact', [0.80, 0.62, 0.42])):
        part = Clip(Input=CleantoGrid(Input=full_cell(block(r, name))))   # merge points on the mirror seams
        part.ClipType = 'Plane'
        part.ClipType.Origin = [0, 0, geo_clip]
        part.ClipType.Normal = [0, 0, 1]
        part.Invert = 1
        surf = ExtractSurface(Input=part)
        d = Show(surf, v)
        d.ColorArrayName = ['POINTS', '']
        d.AmbientColor = d.DiffuseColor = rgb
        d.Specular = 0.25
        d.Ambient = 0.35
        d.Diffuse = 0.75
        d.Interpolation = 'Gouraud'
        outline(part, v, color=(0.1, 0.1, 0.1), width=2.5)
    iso_camera(v, [0.002, -0.002, geo_clip * 0.62], 0.095, (0.55, -0.95, 0.9))
    save(v, 'pele_geometry.png')
    Delete(v)

    # 7: mesh, in the style of a Cubit screenshot (shaded blocks, black edges)
    v = view()
    v.ViewTime = t_end
    for name, rgb in (('matrix_graphite', [0.64, 0.67, 0.76]), ('fuel_compact', [0.86, 0.70, 0.52])):
        d = Show(bottom_slab(full_cell(block(r, name))), v)
        d.SetRepresentationType('Surface With Edges')
        d.ColorArrayName = ['POINTS', '']
        d.AmbientColor = d.DiffuseColor = rgb
        d.EdgeColor = [0.05, 0.05, 0.05]
    iso_camera(v, [0.006, 0.000, SLAB], 0.016, (0.45, -0.85, 1.6))
    save(v, 'pele_mesh_zoomed.png')
    Delete(v)

    v = view()
    v.ViewTime = t_end
    for name, rgb in (('matrix_graphite', [0.64, 0.67, 0.76]), ('fuel_compact', [0.86, 0.70, 0.52])):
        d = Show(bottom_slab(full_cell(block(r, name))), v)
        d.SetRepresentationType('Surface With Edges')
        d.ColorArrayName = ['POINTS', '']
        d.AmbientColor = d.DiffuseColor = rgb
        d.EdgeColor = [0.05, 0.05, 0.05]
    iso_camera(v, WAFER['center'], WAFER['dist'], WAFER['elev'])
    save(v, 'pele_mesh.png')
    Delete(v)


if __name__ == '__main__':
    main()
