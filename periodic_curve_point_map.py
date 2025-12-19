#!/usr/bin/env python3
"""
Periodic unit-cube demo with explicit **point**, **curve**, and **surface** periodic mappings.
- Maps left<->right (x-translation), front<->back (y-translation), bottom<->top (z-translation)
- Applies periodicity on 0D/1D/2D **before** transfinite constraints
- Optionally enforces matching transfinite counts per axis
- Verifies periodic nodes via getPeriodicNodes for surfaces/curves/points

Usage examples:
  python periodic_curve_point_map.py --outfile periodic_verified.msh
  python periodic_curve_point_map.py --transfinite --nx 11 --ny 9 --nz 7 --outfile periodic_verified_tf.msh
"""
import argparse
import gmsh

EPS = 1e-6

# ---------- Utilities for geometry inspection ----------

def get_axis_aligned_faces():
    faces = {}
    for dim, tag in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        if abs(xmin - xmax) < EPS:  # plane x = const
            if abs(xmin) < EPS:
                faces['left'] = tag
            elif abs(xmin - 1.0) < EPS:
                faces['right'] = tag
        elif abs(ymin - ymax) < EPS:  # plane y = const
            if abs(ymin) < EPS:
                faces['front'] = tag
            elif abs(ymin - 1.0) < EPS:
                faces['back'] = tag
        elif abs(zmin - zmax) < EPS:  # plane z = const
            if abs(zmin) < EPS:
                faces['bottom'] = tag
            elif abs(zmin - 1.0) < EPS:
                faces['top'] = tag
    return faces


def entity_bbox(dim, tag):
    return gmsh.model.getBoundingBox(dim, tag)


def on_face_point(pt_tag, face_name):
    x, y, z = gmsh.model.getValue(0, pt_tag, [])
    if face_name == 'left':
        return abs(x - 0.0) < EPS
    if face_name == 'right':
        return abs(x - 1.0) < EPS
    if face_name == 'front':
        return abs(y - 0.0) < EPS
    if face_name == 'back':
        return abs(y - 1.0) < EPS
    if face_name == 'bottom':
        return abs(z - 0.0) < EPS
    if face_name == 'top':
        return abs(z - 1.0) < EPS
    return False


def on_face_curve(c_tag, face_name):
    xmin, ymin, zmin, xmax, ymax, zmax = entity_bbox(1, c_tag)
    if face_name in ('left','right'):
        return abs(xmin - xmax) < EPS and ((face_name=='left' and abs(xmin-0.0)<EPS) or (face_name=='right' and abs(xmin-1.0)<EPS))
    if face_name in ('front','back'):
        return abs(ymin - ymax) < EPS and ((face_name=='front' and abs(ymin-0.0)<EPS) or (face_name=='back' and abs(ymin-1.0)<EPS))
    if face_name in ('bottom','top'):
        return abs(zmin - zmax) < EPS and ((face_name=='bottom' and abs(zmin-0.0)<EPS) or (face_name=='top' and abs(zmin-1.0)<EPS))
    return False


def list_points_on_face(face_name):
    pts = [tag for (d, tag) in gmsh.model.getEntities(0) if on_face_point(tag, face_name)]
    # sort deterministically by in-plane coordinates
    def key(pt):
        x, y, z = gmsh.model.getValue(0, pt, [])
        if face_name in ('left','right'):
            return (y, z)
        if face_name in ('front','back'):
            return (x, z)
        if face_name in ('bottom','top'):
            return (x, y)
        return (x, y, z)
    return sorted(pts, key=key)


def list_curves_on_face(face_name):
    curves = [tag for (d, tag) in gmsh.model.getEntities(1) if on_face_curve(tag, face_name)]
    # sort by mid-point of endpoints in in-plane vars
    def mid_key(c):
        b = gmsh.model.getBoundary([(1, c)], oriented=False, recursive=False)
        p = [t for (d, t) in b if d==0]
        if len(p) >= 2:
            (x1,y1,z1) = gmsh.model.getValue(0, p[0], [])
            (x2,y2,z2) = gmsh.model.getValue(0, p[-1], [])
            xm, ym, zm = (0.5*(x1+x2), 0.5*(y1+y2), 0.5*(z1+z2))
        else:
            xm, ym, zm = gmsh.model.getValue(1, c, [0.5])
        if face_name in ('left','right'):
            return (ym, zm)
        if face_name in ('front','back'):
            return (xm, zm)
        if face_name in ('bottom','top'):
            return (xm, ym)
        return (xm, ym, zm)
    return sorted(curves, key=mid_key)


def affine_translation(tx, ty, tz):
    return [1,0,0,tx, 0,1,0,ty, 0,0,1,tz, 0,0,0,1]

# ---------- Transfinite helpers ----------

def set_transfinite_on_face(face_name, face_tag, nx, ny, nz):
    axis_to_n = {'x': nx, 'y': ny, 'z': nz}
    # set per-curve
    boundary = gmsh.model.getBoundary([(2, face_tag)], oriented=True, recursive=False)
    curves = [t for (d,t) in boundary if d==1]
    for c in curves:
        # deduce main direction of the curve
        b = gmsh.model.getBoundary([(1, c)], oriented=True, recursive=False)
        p = [t for (d,t) in b if d==0]
        (x1,y1,z1) = gmsh.model.getValue(0, p[0], [])
        (x2,y2,z2) = gmsh.model.getValue(0, p[-1], [])
        dx,dy,dz = abs(x2-x1), abs(y2-y1), abs(z2-z1)
        axis = 'x' if dx>=dy and dx>=dz else ('y' if dy>=dx and dy>=dz else 'z')
        gmsh.model.mesh.setTransfiniteCurve(c, axis_to_n[axis], 'Progression', 1.0)
    gmsh.model.mesh.setTransfiniteSurface(face_tag)
    gmsh.model.mesh.setRecombine(2, face_tag)

# ---------- Verification ----------

def verify_periodic_nodes1(dim, target_tags):
    for t in target_tags:
        try:
            mTag, nodes, T = gmsh.model.mesh.getPeriodicNodes(dim, t)
            print(f"[verify] dim={dim} target={t}: master={mTag}, periodic_nodes={len(nodes)}")
        except Exception as e:
            print(f"[verify] dim={dim} target={t}: getPeriodicNodes failed: {e}")


def verify_periodic_nodes(dim, target_tags):
    for t in target_tags:
        try:
            res = gmsh.model.mesh.getPeriodicNodes(dim, t)
            # API differences across versions: res can be 3-tuple or 4-tuple
            if isinstance(res, (tuple, list)):
                if len(res) == 3:
                    masterTag, nodes, T = res
                    print(f"[verify] dim={dim} target={t}: master={masterTag}, periodic_nodes={len(nodes)}")
                elif len(res) == 4:
                    masterTag, nodes, nodeMap, T = res
                    print(f"[verify] dim={dim} target={t}: master={masterTag}, periodic_nodes={len(nodes)}, map_len={len(nodeMap)}")
                else:
                    print(f"[verify] dim={dim} target={t}: unexpected len(res)={len(res)}")
            else:
                print(f"[verify] dim={dim} target={t}: unexpected return type {type(res)}")
        except Exception as e:
            print(f"[verify] dim={dim} target={t}: getPeriodicNodes failed: {e}")

# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--transfinite', action='store_true', help='Apply transfinite on faces with axis-aware counts')
    ap.add_argument('--nx', type=int, default=11, help='Points on x-directed edges')
    ap.add_argument('--ny', type=int, default=9, help='Points on y-directed edges')
    ap.add_argument('--nz', type=int, default=7, help='Points on z-directed edges')
    ap.add_argument('--outfile', type=str, default='periodic_verified.msh', help='Output mesh file')
    args = ap.parse_args()

    gmsh.initialize()
    try:
        gmsh.model.add('periodic_curve_point_map')
        gmsh.model.occ.addBox(0,0,0, 1,1,1)
        gmsh.model.occ.synchronize()

        faces = get_axis_aligned_faces()
        # Prepare periodic pairs and translations
        pairs = [
            ('right','left',  (1,0,0)),
            ('back','front',  (0,1,0)),
            ('top','bottom',  (0,0,1)),
        ]

        # 1) Periodicity on POINTS and CURVES (before any transfinite)
        all_points = [t for (d,t) in gmsh.model.getEntities(0)]
        all_curves = [t for (d,t) in gmsh.model.getEntities(1)]

        for tgt_name, src_name, v in pairs:
            # Points
            tgt_pts = list_points_on_face(tgt_name)
            src_pts = list_points_on_face(src_name)
            if len(tgt_pts) != len(src_pts):
                print(f"[warn] point count mismatch {tgt_name}({len(tgt_pts)}) vs {src_name}({len(src_pts)})")
            if tgt_pts and src_pts:
                gmsh.model.mesh.setPeriodic(0, tgt_pts, src_pts, affine_translation(*v))

            # Curves
            tgt_cur = list_curves_on_face(tgt_name)
            src_cur = list_curves_on_face(src_name)
            if len(tgt_cur) != len(src_cur):
                print(f"[warn] curve count mismatch {tgt_name}({len(tgt_cur)}) vs {src_name}({len(src_cur)})")
            if tgt_cur and src_cur:
                gmsh.model.mesh.setPeriodic(1, tgt_cur, src_cur, affine_translation(*v))

        # 2) Periodicity on SURFACES
        for tgt_name, src_name, v in pairs:
            gmsh.model.mesh.setPeriodic(2, [faces[tgt_name]], [faces[src_name]], affine_translation(*v))

        # 3) Now (optionally) apply transfinite on each face with axis-aware counts
        if args.transfinite:
            set_transfinite_on_face('left',   faces['left'],   args.nx, args.ny, args.nz)
            set_transfinite_on_face('right',  faces['right'],  args.nx, args.ny, args.nz)
            set_transfinite_on_face('front',  faces['front'],  args.nx, args.ny, args.nz)
            set_transfinite_on_face('back',   faces['back'],   args.nx, args.ny, args.nz)
            set_transfinite_on_face('bottom', faces['bottom'], args.nx, args.ny, args.nz)
            set_transfinite_on_face('top',    faces['top'],    args.nx, args.ny, args.nz)

        # Size hints (optional, just to visualize non-uniformity)
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.12)
        gmsh.model.mesh.setSize([(0,1)], 0.02)

        # Generate and write
        gmsh.model.mesh.generate(3)
        gmsh.write(args.outfile)
        print(f"Wrote mesh to {args.outfile}")

        # 4) Verify periodic nodes on faces/curves/points
        verify_periodic_nodes(2, [faces['right'], faces['back'], faces['top']])
        # For curves and points, just check a few representatives (all targets on right/back/top)
        verify_periodic_nodes(1, list_curves_on_face('right'))
        verify_periodic_nodes(0, list_points_on_face('right'))

    finally:
        gmsh.finalize()

if __name__ == '__main__':
    main()
