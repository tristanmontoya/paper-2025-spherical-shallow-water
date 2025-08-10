# === Axis-hugging geo labels (uses *display* bounds) — BLACK ===
# Longitudes: 150°W .. 150°E (step 30) along X at y=ymin
# Latitudes :  60°S  .. 60°N  (step 30) along Y at x=xmin

from paraview.simple import *
from math import sin, cos, radians

# ticks
LON_MIN, LON_MAX, LON_STEP = -150, 150, 30
LAT_MIN, LAT_MAX, LAT_STEP =  -60,  60, 30

# style
FONT_SIZE = 18
COLOR = [0.0, 0.0, 0.0]     # black
PREFIX = "Geo3D"

# nudge labels away from axes (as fraction of axis length)
XAXIS_Y_OFFSET = 0.02
YAXIS_X_OFFSET = 0.02

rv = GetActiveViewOrCreate('RenderView')
src = GetActiveSource()
if not src:
    raise RuntimeError("Select your dataset in the Pipeline Browser first.")

# ensure it's shown and pipeline updated
rep = GetDisplayProperties(src, view=rv)
Show(src, rv)
Render()

def cleanup():
    for (name, _id), p in list(GetSources().items()):
        if str(name).startswith(PREFIX + "_"):
            try: Delete(p)
            except: pass

def fmt_lon(d):
    return f"{abs(d)}\N{DEGREE SIGN} W" if d < 0 else (f"{d}\N{DEGREE SIGN} E" if d > 0 else "0\N{DEGREE SIGN}")

def fmt_lat(d):
    return f"{abs(d)}\N{DEGREE SIGN} S" if d < 0 else (f"{d}\N{DEGREE SIGN} N" if d > 0 else "0\N{DEGREE SIGN}")

def rotM(ang):
    # Euler XYZ (degrees) -> 3x3
    ax, ay, az = [radians(a) for a in ang]
    cx, sx = cos(ax), sin(ax)
    cy, sy = cos(ay), sin(ay)
    cz, sz = cos(az), sin(az)
    Rx = ((1,0,0),(0,cx,-sx),(0,sx,cx))
    Ry = ((cy,0,sy),(0,1,0),(-sy,0,cy))
    Rz = ((cz,-sz,0),(sz,cz,0),(0,0,1))
    # R = Rz * Ry * Rx
    def matmul(A,B):
        return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)) for i in range(3))
    return matmul(Rz, matmul(Ry, Rx))

def xform_pt(p, origin, scale, orient, pos):
    # (p-origin)->scale->rotate->( + origin + pos )
    px,py,pz = p[0]-origin[0], p[1]-origin[1], p[2]-origin[2]
    px,py,pz = px*scale[0], py*scale[1], pz*scale[2]
    R = rotM(orient)
    rx = R[0][0]*px + R[0][1]*py + R[0][2]*pz
    ry = R[1][0]*px + R[1][1]*py + R[1][2]*pz
    rz = R[2][0]*px + R[2][1]*py + R[2][2]*pz
    return (rx + origin[0] + pos[0],
            ry + origin[1] + pos[1],
            rz + origin[2] + pos[2])

def display_bounds(r):
    # get rendered (represented) data bounds
    try:
        di = r.GetRepresentedDataInformation()
    except Exception:
        di = None
    if di is None:
        di = src.GetDataInformation()
    b = list(di.GetBounds())
    xmin,xmax,ymin,ymax,zmin,zmax = b
    # apply Display transform
    def getv(obj, name, default):
        try: return list(getattr(obj, name))
        except Exception: return list(default)
    pos    = getv(r, 'Position',    [0,0,0])
    scale  = getv(r, 'Scale',       [1,1,1])
    orient = getv(r, 'Orientation', [0,0,0])
    origin = getv(r, 'Origin',      [0,0,0])

    corners = [(xmin,ymin,zmin),(xmin,ymin,zmax),(xmin,ymax,zmin),(xmin,ymax,zmax),
               (xmax,ymin,zmin),(xmax,ymin,zmax),(xmax,ymax,zmin),(xmax,ymax,zmax)]
    tx,ty,tz = zip(*[xform_pt(c, origin, scale, orient, pos) for c in corners])
    return (min(tx), max(tx), min(ty), max(ty), min(tz), max(tz))

def lerp(v, a0, a1, b0, b1):
    if a1 == a0: return 0.5*(b0+b1)
    t = (v - a0) / float(a1 - a0)
    t = max(0.0, min(1.0, t))
    return b0 + t*(b1 - b0)

def ensure_billboard(disp):
    try:
        disp.TextPropMode = 'Billboard 3D Text'
    except Exception:
        pass

def add_label(name, text, pos_xyz, h_just='Centered', v_just='Centered'):
    t = Text(registrationName=name, Text=text)
    d = Show(t, rv)
    ensure_billboard(d)
    try: d.BillboardPosition = pos_xyz
    except Exception:
        # Fallback to 3D geometry if billboard is not available in your build
        vg = VectorText(Text=text)
        tr = Transform(Input=vg)
        tr.Transform = 'Transform'
        tr.Transform.Translate = pos_xyz
        t = tr
        d = Show(t, rv)
    # style (guard per-version differences)
    for (prop,val) in [('Justification',h_just), ('VerticalJustification',v_just),
                       ('FontSize',FONT_SIZE), ('Color',COLOR), ('TextColor',COLOR)]:
        try: setattr(d, prop, val)
        except Exception: pass
    return t

# ---------- build ----------
cleanup()

xmin,xmax,ymin,ymax,zmin,zmax = display_bounds(rep)
Lx = xmax - xmin
Ly = ymax - ymin

# Bottom: longitudes along X at y=ymin - offset
y_bottom = ymin - XAXIS_Y_OFFSET * Ly
z_plane  = zmin
for lon in range(LON_MIN, LON_MAX+1, LON_STEP):
    x = lerp(lon, LON_MIN, LON_MAX, xmin, xmax)
    add_label(f"{PREFIX}_lon_{lon:+03d}", fmt_lon(lon), [x, y_bottom, z_plane],
              h_just='Centered', v_just='Top')

# Left: latitudes along Y at x=xmin - offset
x_left = xmin - YAXIS_X_OFFSET * Lx
for lat in range(LAT_MIN, LAT_MAX+1, LAT_STEP):
    y = lerp(lat, LAT_MIN, LAT_MAX, ymin, ymax)
    add_label(f"{PREFIX}_lat_{lat:+03d}", fmt_lat(lat), [x_left, y, z_plane],
              h_just='Right', v_just='Centered')

Render()
# === end ===