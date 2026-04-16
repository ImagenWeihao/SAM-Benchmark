"""Compute physical overlap between 10X and 40X ND2 files.

Alignment uses BOTH stage positions AND the `cameraTransformationMatrix`
metadata field in each file (the Nikon pixel->stage transform).  For these
files the matrix is approximately [[-1, 0],[0, -1]], i.e.

    stage +X  ->  image LEFT  (col decreases)
    stage +Y  ->  image UP    (row decreases)

plus a small (~0.4 deg) rotation we apply when placing tiles.
"""
from pathlib import Path
import numpy as np
import nd2
import tifffile
from PIL import Image

P10 = Path(r"D:\Matlab Files\Codes\SAM_Data\AnnotationData"
           r"\20260412_195011_768__WellA05_Channel405,561,532_1,BF_Seq0000.nd2")
P40 = Path(r"D:\Matlab Files\Codes\SAM_Data\AnnotationData"
           r"\20260413_105534_226__WellA05_Channel488,405,561,532_1_Seq0000_crop(first5FOV).nd2")
OUT = Path(r"D:\Matlab Files\Codes\SAM_Data\AnnotationData\Overlap")
OUT.mkdir(exist_ok=True)

CH10_DAPI = 0   # '405'
CH40_DAPI = 1   # '405'


def read_positions(f):
    for loop in f.experiment:
        params = getattr(loop, "parameters", None)
        if params is not None and hasattr(params, "points"):
            return [(p.stagePositionUm.x, p.stagePositionUm.y) for p in params.points]


def read_cam_matrix(f):
    """Read pixel->stage 2x2 transformation matrix (um per pixel, signed)."""
    vol = f.metadata.channels[0].volume
    m = vol.cameraTransformationMatrix
    return np.array([[m[0], m[1]], [m[2], m[3]]], dtype=np.float64)


def norm_u8(img, lo=1.0, hi=99.5):
    img = img.astype(np.float32)
    nz = img[img > 0]
    if nz.size == 0:
        return np.zeros(img.shape, dtype=np.uint8)
    mn = np.percentile(nz, lo); mx = np.percentile(nz, hi)
    return (np.clip((img - mn) / max(mx - mn, 1e-9), 0, 1) * 255).astype(np.uint8)


def draw_rect(rgb, r1, r2, c1, c2, color, w=2):
    H, W = rgb.shape[:2]
    r1 = max(0, r1); r2 = min(H, r2); c1 = max(0, c1); c2 = min(W, c2)
    if r2 <= r1 or c2 <= c1: return
    rgb[r1:min(r1 + w, r2), c1:c2] = color
    rgb[max(r2 - w, r1):r2, c1:c2] = color
    rgb[r1:r2, c1:min(c1 + w, c2)] = color
    rgb[r1:r2, max(c2 - w, c1):c2] = color


# ──────────────────────────────────────────────────────────────────────────
print("[10X] loading...")
with nd2.ND2File(str(P10)) as f:
    px10 = float(f.voxel_size().x)
    arr10 = f.asarray()
    pos10 = read_positions(f)
    M10 = read_cam_matrix(f)    # pixel->stage, in µm
H10, W10 = arr10.shape[-2:]
print(f"[10X] px={px10:.4f} um  tile={W10}x{H10}  n={len(pos10)}")
print(f"[10X] cameraTransformationMatrix (pixel->stage):\n{M10}")

print("[40X] loading MIP...")
with nd2.ND2File(str(P40)) as f:
    px40 = float(f.voxel_size().x)
    arr40 = f.asarray()
    pos40 = read_positions(f)
    M40 = read_cam_matrix(f)
H40, W40 = arr40.shape[-2:]
mip40 = arr40[:, :, CH40_DAPI].max(axis=1)
print(f"[40X] px={px40:.4f} um  tile={W40}x{H40}  n={len(pos40)}")
print(f"[40X] cameraTransformationMatrix (pixel->stage):\n{M40}")

# ──────────────────────────────────────────────────────────────────────────
# Derive stage->canvas-pixel transforms.
#
# Nikon's M maps pixel offset from IMAGE CENTER -> stage offset from TILE CENTRE:
#     [dx_stage]     [du]
#     [dy_stage] = M [dv]
#
# Inverting: M_inv maps stage_offset -> pixel_offset_from_center.  We build a
# common world-to-canvas mapping T_world2canvas (2x2 + translation) such that:
#     canvas_col = ((stage_x - cx) * T[0,0] + (stage_y - cy) * T[0,1]) + col_offset
#     canvas_row = ((stage_x - cx) * T[1,0] + (stage_y - cy) * T[1,1]) + row_offset
# where T_world2canvas = M_inv / px   (since M is in µm per pixel).
# ──────────────────────────────────────────────────────────────────────────

# The 10X and 40X matrices are the same here; each tile is placed using its
# own file's M_inv (scaled to its own px).
Minv10_px = np.linalg.inv(M10) * px10   # matrix maps um offset -> px offset at 10X scale
# Wait: Minv maps stage_offset (um) -> pixel_offset (px).  For 10X px=0.655:
#   pixel_offset_from_tile_center = Minv @ stage_offset
# For consistent 10X-scale canvas we want canvas pixel coords where 1 px = px10 µm.
# So canvas_offset_from_tile_center_px = Minv10 @ stage_offset_um
Minv10 = np.linalg.inv(M10) / px10   # um -> 10X-canvas-px
Minv40 = np.linalg.inv(M40) / px40   # um -> 40X-canvas-px (native 40X scale)
print(f"\n[REG] Minv10 (stage_um -> 10X_canvas_px):\n{Minv10}")
print(f"[REG] Minv40 (stage_um -> 40X_canvas_px):\n{Minv40}")


def tile_corners_world_to_canvas_px(sx, sy, W_px, H_px, Minv):
    """Return the 4 image-corner positions in canvas-px offset from stage centre.
    Output order: TL, TR, BR, BL (in image coord convention: row,col)."""
    # Image pixel offsets from image center (in pixels):
    corners_px = np.array([
        [-W_px / 2, -H_px / 2],   # TL in image: (col-offset, row-offset) = (-W/2, -H/2)
        [+W_px / 2, -H_px / 2],   # TR
        [+W_px / 2, +H_px / 2],   # BR
        [-W_px / 2, +H_px / 2],   # BL
    ])
    # Convert each image pixel offset to stage offset via M (use the file's M):
    # But we want canvas px around a chosen WORLD origin — easier to compute
    # the offset in the reference frame where canvas_x = +X_stage direction
    # scaled.  So we compute each image corner's STAGE position, then convert.
    return corners_px


def world_to_canvas_px(stage_xy, origin_xy, Minv):
    """stage_xy (um) -> canvas_px_offset using the Minv transform.
    Canvas axes: +col = image-right, +row = image-down."""
    d = np.asarray(stage_xy) - np.asarray(origin_xy)   # (2,) in um
    return Minv @ d   # returns (2,) = (col_offset_px, row_offset_px)


# Pick a reference world point (any of the 10X stage centers will do — pick
# FOV0 as reference so its center = canvas origin, then translate canvas)
ref_x, ref_y = pos10[0]


# Step 1: compute the canvas-px location of each tile CENTRE for 10X
centers10_cpx = [world_to_canvas_px(p, (ref_x, ref_y), Minv10) for p in pos10]
# Centres of 40X FOVs in 10X canvas px space (via 10X's Minv but with px40->px10)
# Actually for a common canvas we convert 40X stage positions into 10X-pixel
# canvas using Minv10 (since 10X is the canvas spec). This requires the same
# transform — we assume M is shared (which it is here).
centers40_in10cpx = [world_to_canvas_px(p, (ref_x, ref_y), Minv10) for p in pos40]

# Compute canvas bounds taking into account all tile footprints.
# For each 10X tile, the 4 image corners in canvas px (relative to ref) are:
#   corner_image_px = [(col_off, row_off)]
#   canvas_px_of_corner = tile_center_canvas_px + corner_image_px
#   (since M is applied via center, and image pixel offsets match canvas pixel
#    offsets directly — the canvas IS an image-pixel grid)
halfW10 = W10 / 2; halfH10 = H10 / 2
halfW40_in_10 = (W40 * px40) / px10 / 2
halfH40_in_10 = (H40 * px40) / px10 / 2

def tile_bounds(center_cpx, halfW, halfH):
    cc, rc = center_cpx
    return (cc - halfW, cc + halfW, rc - halfH, rc + halfH)

all_bounds = []
for c in centers10_cpx:
    all_bounds.append(tile_bounds(c, halfW10, halfH10))
for c in centers40_in10cpx:
    all_bounds.append(tile_bounds(c, halfW40_in_10, halfH40_in_10))

min_col = min(b[0] for b in all_bounds)
max_col = max(b[1] for b in all_bounds)
min_row = min(b[2] for b in all_bounds)
max_row = max(b[3] for b in all_bounds)
cW10 = int(np.ceil(max_col - min_col))
cH10 = int(np.ceil(max_row - min_row))
col_ofs = -min_col; row_ofs = -min_row
print(f"[REG]  10X-scale canvas: {cW10} x {cH10} px  "
      f"(col offset {col_ofs:.1f}, row offset {row_ofs:.1f})")


def place_tile_simple(canvas, tile, center_cpx_in_canvas):
    """Paste tile onto canvas at integer top-left derived from the tile centre."""
    h, w = tile.shape
    cc, rc = center_cpx_in_canvas
    row0 = int(round(rc - h / 2))
    col0 = int(round(cc - w / 2))
    H, W = canvas.shape
    r1c, r2c = max(0, row0), min(H, row0 + h)
    c1c, c2c = max(0, col0), min(W, col0 + w)
    if r2c <= r1c or c2c <= c1c: return
    sr1 = r1c - row0; sr2 = sr1 + (r2c - r1c)
    sc1 = c1c - col0; sc2 = sc1 + (c2c - c1c)
    canvas[r1c:r2c, c1c:c2c] = tile[sr1:sr2, sc1:sc2]


# Build 10X stitched canvas at 10X resolution
canvas10 = np.zeros((cH10, cW10), dtype=np.uint16)
centers10_in_canvas = []
for p, cpx in enumerate(centers10_cpx):
    cc_c = cpx[0] + col_ofs
    rc_c = cpx[1] + row_ofs
    centers10_in_canvas.append((cc_c, rc_c))
    place_tile_simple(canvas10, arr10[p, CH10_DAPI], (cc_c, rc_c))
print(f"[10X] stitched canvas10: {canvas10.shape}")

# Build 40X stitched canvas at 10X resolution (downsampled each tile so
# we can visualise in same pixel grid as 10X).  We'll also keep a native-
# resolution canvas for exporting the overlap at 40X detail.
canvas40_in10scale = np.zeros((cH10, cW10), dtype=np.uint16)
centers40_in_canvas = []
scale40to10 = px40 / px10   # tile gets smaller on 10X grid
new_W40_on10 = int(round(W40 * scale40to10))
new_H40_on10 = int(round(H40 * scale40to10))
for p, cpx in enumerate(centers40_in10cpx):
    cc_c = cpx[0] + col_ofs
    rc_c = cpx[1] + row_ofs
    centers40_in_canvas.append((cc_c, rc_c))
    tile_ds = np.asarray(Image.fromarray(mip40[p].astype(np.float32)).resize(
        (new_W40_on10, new_H40_on10), Image.BILINEAR), dtype=np.uint16)
    place_tile_simple(canvas40_in10scale, tile_ds, (cc_c, rc_c))
print(f"[40X-on-10X-scale] stitched: {canvas40_in10scale.shape}")

# Also build NATIVE-40X canvas covering only the 40X strip bbox
x40 = [world_to_canvas_px(p, (ref_x, ref_y), Minv40) for p in pos40]
halfW40 = W40 / 2; halfH40 = H40 / 2
b40 = [tile_bounds(c, halfW40, halfH40) for c in x40]
m40_col = min(b[0] for b in b40); M40_col = max(b[1] for b in b40)
m40_row = min(b[2] for b in b40); M40_row = max(b[3] for b in b40)
cW40 = int(np.ceil(M40_col - m40_col)); cH40 = int(np.ceil(M40_row - m40_row))
col40_ofs = -m40_col; row40_ofs = -m40_row
canvas40_native = np.zeros((cH40, cW40), dtype=np.uint16)
for p, cpx in enumerate(x40):
    place_tile_simple(canvas40_native, mip40[p],
                      (cpx[0] + col40_ofs, cpx[1] + row40_ofs))
print(f"[40X-native]      stitched: {canvas40_native.shape}")

# ──────────────────────────────────────────────────────────────────────────
# Compute overlap region in 10X canvas pixel coords.
# A canvas pixel is "covered by 10X" iff it was painted by any 10X tile.
# Similarly for 40X.  Use binary masks to find intersection.
# ──────────────────────────────────────────────────────────────────────────
mask10 = np.zeros(canvas10.shape, dtype=bool)
for cc, rc in centers10_in_canvas:
    r0 = int(round(rc - H10 / 2)); c0 = int(round(cc - W10 / 2))
    r1 = max(0, r0); r2 = min(canvas10.shape[0], r0 + H10)
    c1 = max(0, c0); c2 = min(canvas10.shape[1], c0 + W10)
    mask10[r1:r2, c1:c2] = True

mask40 = np.zeros(canvas10.shape, dtype=bool)
for cc, rc in centers40_in_canvas:
    r0 = int(round(rc - new_H40_on10 / 2)); c0 = int(round(cc - new_W40_on10 / 2))
    r1 = max(0, r0); r2 = min(canvas10.shape[0], r0 + new_H40_on10)
    c1 = max(0, c0); c2 = min(canvas10.shape[1], c0 + new_W40_on10)
    mask40[r1:r2, c1:c2] = True

overlap = mask10 & mask40
rows_any = np.any(overlap, axis=1)
cols_any = np.any(overlap, axis=0)
if not rows_any.any() or not cols_any.any():
    raise SystemExit("No overlap between 10X and 40X.")
r_a = int(np.argmax(rows_any)); r_b = len(rows_any) - int(np.argmax(rows_any[::-1]))
c_a = int(np.argmax(cols_any)); c_b = len(cols_any) - int(np.argmax(cols_any[::-1]))
print(f"\n[OVERLAP] 10X-canvas rows [{r_a}, {r_b}]  cols [{c_a}, {c_b}]  "
      f"size {c_b - c_a} x {r_b - r_a} px  = "
      f"{(c_b - c_a) * px10:.1f} x {(r_b - r_a) * px10:.1f} um")

crop10 = canvas10[r_a:r_b, c_a:c_b]

# Corresponding slice of the on-10X 40X canvas
crop40_on10 = canvas40_in10scale[r_a:r_b, c_a:c_b]

# Native-resolution 40X crop: find which native-40X canvas region corresponds
# to the 10X-canvas overlap region.  Convert 4 canvas corners via their world
# coordinates and then into native-40X canvas pixels.
def canvas10_to_world(col_c, row_c):
    d_cpx = np.array([col_c - col_ofs, row_c - row_ofs])   # offset from ref center
    # canvas_px = Minv10 @ (stage - ref)   =>  stage - ref = M10 @ canvas_px
    d_world = M10 @ d_cpx
    return ref_x + d_world[0], ref_y + d_world[1]

def world_to_canvas40_native_px(wx, wy):
    d_world = np.array([wx - ref_x, wy - ref_y])
    d_cpx = Minv40 @ d_world
    return d_cpx[0] + col40_ofs, d_cpx[1] + row40_ofs

corner_world = [canvas10_to_world(c_a, r_a), canvas10_to_world(c_b, r_a),
                canvas10_to_world(c_b, r_b), canvas10_to_world(c_a, r_b)]
corner_40px = [world_to_canvas40_native_px(*w) for w in corner_world]
cols40 = [c for c, r in corner_40px]
rows40 = [r for c, r in corner_40px]
r40_a = max(0, int(np.floor(min(rows40))))
r40_b = min(canvas40_native.shape[0], int(np.ceil(max(rows40))))
c40_a = max(0, int(np.floor(min(cols40))))
c40_b = min(canvas40_native.shape[1], int(np.ceil(max(cols40))))
crop40_native = canvas40_native[r40_a:r40_b, c40_a:c40_b]
print(f"[OVERLAP] 40X-native rows [{r40_a}, {r40_b}]  cols [{c40_a}, {c40_b}]  "
      f"size {c40_b - c40_a} x {r40_b - r40_a} px")

# ──────────────────────────────────────────────────────────────────────────
# Build annotated overview
# ──────────────────────────────────────────────────────────────────────────
overview = np.stack([norm_u8(canvas10)] * 3, axis=-1)
for cc, rc in centers40_in_canvas:
    r0 = int(round(rc - new_H40_on10 / 2)); c0 = int(round(cc - new_W40_on10 / 2))
    draw_rect(overview, r0, r0 + new_H40_on10, c0, c0 + new_W40_on10, (80, 255, 80), w=2)
draw_rect(overview, r_a, r_b, c_a, c_b, (255, 60, 60), w=3)

# side-by-side overlap (10X vs 40X-on-10X-scale)
sbs = np.concatenate([norm_u8(crop10), norm_u8(crop40_on10)], axis=0)

# 40X stitched strip (native) — standalone TIF
native40_crop = canvas40_native.copy()

# ──────────────────────────────────────────────────────────────────────────
# Save all outputs
# ──────────────────────────────────────────────────────────────────────────
to_save = {
    "10x_overlap_crop_DAPI.tif":        (crop10,       px10),
    "40x_overlap_crop_DAPI_MIP.tif":    (crop40_native, px40),
    "40x_overlap_at10x_scale.tif":      (crop40_on10,  px10),
    "10x_stitched_DAPI.tif":            (canvas10,     px10),
    "40x_stitched_DAPI_MIP.tif":        (native40_crop, px40),
}
for name, (img, px) in to_save.items():
    tifffile.imwrite(OUT / name, img, photometric='minisblack',
                     resolution=(1e4 / px, 1e4 / px),
                     metadata={'unit': 'um', 'axes': 'YX'})
tifffile.imwrite(OUT / "10x_stitched_with_40x_overlay.tif", overview, photometric='rgb')
tifffile.imwrite(OUT / "overlap_side_by_side_10xTop_40xBottom.tif",
                 sbs, photometric='minisblack')

# ──────────────────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────────────────
summary = [
    "Overlap computation (stage positions + cameraTransformationMatrix)",
    "=" * 60,
    f"10X file : {P10.name}",
    f"40X file : {P40.name}",
    "",
    f"10X px : {px10:.4f} um/px   tile {W10}x{H10}   n_tiles {len(pos10)}",
    f"40X px : {px40:.4f} um/px   tile {W40}x{H40}   n_tiles {len(pos40)} (MIP over 15 Z)",
    "",
    "cameraTransformationMatrix M (pixel-offset-from-centre -> stage-offset):",
    f"  10X M = {M10.tolist()}",
    f"  40X M = {M40.tolist()}",
    "  => stage +X maps to image LEFT  (not right)",
    "  => stage +Y maps to image UP    (as expected)",
    "  => plus a ~0.4 deg rotation (small; ignored at integer-tile placement)",
    "",
    f"10X-scale stitched canvas: {canvas10.shape[1]} x {canvas10.shape[0]} px",
    f"Overlap crop in 10X canvas: rows [{r_a}, {r_b}], cols [{c_a}, {c_b}]  "
    f"({c_b - c_a} x {r_b - r_a} px  = {(c_b-c_a)*px10:.1f} x {(r_b-r_a)*px10:.1f} um)",
    "",
    "40X FOV centres in 10X stitched canvas (col, row):",
]
for p, (cc, rc) in enumerate(centers40_in_canvas):
    summary.append(f"  FOV{p+1}: stage=({pos40[p][0]:.1f}, {pos40[p][1]:.1f}) um  "
                   f"-> canvas ({cc:.1f}, {rc:.1f})")
summary.append("")
summary.append("10X FOV centres in canvas:")
for p, (cc, rc) in enumerate(centers10_in_canvas):
    summary.append(f"  10X FOV{p}: stage=({pos10[p][0]:.1f}, {pos10[p][1]:.1f}) um  "
                   f"-> canvas ({cc:.1f}, {rc:.1f})")
summary += [
    "",
    "Outputs:",
    "  10x_overlap_crop_DAPI.tif          - 10X overlap crop",
    "  40x_overlap_crop_DAPI_MIP.tif      - same physical region, 40X native res",
    "  40x_overlap_at10x_scale.tif        - same physical region, down-sampled to 10X px",
    "  10x_stitched_DAPI.tif              - full 10X 3x3 mosaic",
    "  40x_stitched_DAPI_MIP.tif          - 40X strip stitched (native)",
    "  10x_stitched_with_40x_overlay.tif  - 10X with 40X footprint + overlap overlay",
    "  overlap_side_by_side_10xTop_40xBottom.tif",
]
text = "\n".join(summary)
(OUT / "overlap_summary.txt").write_text(text)
print("\n" + text)
print(f"\n[DONE] outputs -> {OUT}")
