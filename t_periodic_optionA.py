#!/usr/bin/env python3
import argparse
import gmsh

EPS = 1e-5

# ----------------------- Geometry helpers -----------------------

def find_axis_aligned_faces():
    """Return dict with keys 'left','right','front','back','bottom','top' -> surface tags.
    Uses bounding boxes to identify faces of the unit box created at (0,0,0)->(1,1,1).
    """
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


def get_surface_edges(surf_tag):
    b = gmsh.model.getBoundary([(2, surf_tag)], oriented=True, recursive=False)
    return [t for d, t in b if d == 1]


def edge_axis(curve_tag):
    pts = gmsh.model.getBoundary([(1, curve_tag)], oriented=True, recursive=False)
    end_tags = [t for d, t in pts if d == 0]
    xyz = [gmsh.model.getValue(0, p, []) for p in end_tags]
    dx = abs(xyz[0][0] - xyz[1][0])
    dy = abs(xyz[0][1] - xyz[1][1])
    dz = abs(xyz[0][2] - xyz[1][2])
    if dx >= dy and dx >= dz:
        return 'x'
    elif dy >= dx and dy >= dz:
        return 'y'
    else:
        return 'z'

# ----------------------- Periodic helpers -----------------------

def affine_translation(tx, ty, tz):
    return [1,0,0,tx,
            0,1,0,ty,
            0,0,1,tz,
            0,0,0,1]


def list_points_on_named_face(face_name, face_tag):
    """Return *sorted* list of point tags classified on a face, sorted by in-plane coords
    to make the target/source lists correspond when calling setPeriodic(0,...).
    """
    # collect all boundary points of the face
    boundary_curves = gmsh.model.getBoundary([(2, face_tag)], oriented=False, recursive=False)
    pts = []
    for (d, ctag) in boundary_curves:
        if d != 1:
            continue
        curve_pts = gmsh.model.getBoundary([(1, ctag)], oriented=False, recursive=False)
        pts.extend([t for (dd, t) in curve_pts if dd == 0])
    unique_pts = sorted(set(pts))

    def key(pt):
        x, y, z = gmsh.model.getValue(0, pt, [])
        if face_name in ('left','right'):
            return (y, z)
        if face_name in ('front','back'):
            return (x, z)
        if face_name in ('bottom','top'):
            return (x, y)
        return (x, y, z)

    return sorted(unique_pts, key=key)


def list_curves_on_named_face(face_name, face_tag):
    """Return *sorted* list of curve tags lying on a given face.
    Sorting by in-plane mid-point for stable target/source ordering.
    """
    curves = [t for (d, t) in gmsh.model.getBoundary([(2, face_tag)], oriented=False, recursive=False) if d == 1]

    def mid_key(c):
        b = gmsh.model.getBoundary([(1, c)], oriented=False, recursive=False)
        p = [t for (d, t) in b if d == 0]
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


def set_periodic_surface(target_tag, source_tag, translation_vec):
    tx, ty, tz = translation_vec
    T = affine_translation(tx, ty, tz)
    gmsh.model.mesh.setPeriodic(2, [target_tag], [source_tag], T)

# ----------------------- Transfinite helpers -----------------------

def set_transfinite_rect_surface(surf_tag, nx, ny):
    edges = get_surface_edges(surf_tag)
    x_edges = []
    y_edges = []
    for c in edges:
        ax = edge_axis(c)
        if ax == 'x':
            x_edges.append(c)
        elif ax == 'y':
            y_edges.append(c)
    for c in x_edges:
        gmsh.model.mesh.setTransfiniteCurve(c, nx, 'Progression', 1.0)
    for c in y_edges:
        gmsh.model.mesh.setTransfiniteCurve(c, ny, 'Progression', 1.0)

    
    gmsh.model.mesh.setTransfiniteSurface(surf_tag)
    gmsh.model.mesh.setRecombine(2, surf_tag)

# ----------------------- Verification -----------------------

def verify_periodic_nodes(dim, target_tags):
    for t in target_tags:
        try:
            res = gmsh.model.mesh.getPeriodicNodes(dim, t)
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

# ----------------------- Main -----------------------

def main():
    parser = argparse.ArgumentParser(description='Periodic + Transfinite demo (Option A: map points & curves before transfinite)')
    parser.add_argument('--transfinite', action='store_true', help='Enable transfinite quad faces with matching counts')
    parser.add_argument('--nx', type=int, default=9, help='#points along x-directed edges of faces (in-plane x)')
    parser.add_argument('--ny', type=int, default=7, help='#points along y-directed edges of faces (in-plane y)')
    parser.add_argument('--outfile', type=str, default='t_periodic_optionA.msh', help='Output mesh file')
    args = parser.parse_args()

    gmsh.initialize()
    try:
        gmsh.model.add('t_periodic_optionA')
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
        gmsh.model.occ.synchronize()

        # Size hints (optional)
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.12)
        gmsh.model.mesh.setSize([(0, 1)], 0.02)

        faces = find_axis_aligned_faces()

        # ---- Periodicity on surfaces, curves, and points (in this order) ----
        # right <- left (x+1)
        set_periodic_surface(faces['right'], faces['left'], (1,0,0))
        gmsh.model.mesh.setPeriodic(1,
            list_curves_on_named_face('right', faces['right']),
            list_curves_on_named_face('left',  faces['left']),
            affine_translation(1,0,0))
        gmsh.model.mesh.setPeriodic(0,
            list_points_on_named_face('right', faces['right']),
            list_points_on_named_face('left',  faces['left']),
            affine_translation(1,0,0))

        # back <- front (y+1)
        set_periodic_surface(faces['back'], faces['front'], (0,1,0))
        gmsh.model.mesh.setPeriodic(1,
            list_curves_on_named_face('back',  faces['back']),
            list_curves_on_named_face('front', faces['front']),
            affine_translation(0,1,0))
        gmsh.model.mesh.setPeriodic(0,
            list_points_on_named_face('back',  faces['back']),
            list_points_on_named_face('front', faces['front']),
            affine_translation(0,1,0))

        # top <- bottom (z+1)
        set_periodic_surface(faces['top'], faces['bottom'], (0,0,1))
        gmsh.model.mesh.setPeriodic(1,
            list_curves_on_named_face('top',    faces['top']),
            list_curves_on_named_face('bottom', faces['bottom']),
            affine_translation(0,0,1))
        gmsh.model.mesh.setPeriodic(0,
            list_points_on_named_face('top',    faces['top']),
            list_points_on_named_face('bottom', faces['bottom']),
            affine_translation(0,0,1))

        # ---- Now (optionally) apply transfinite constraints ----
        if args.transfinite:
            set_transfinite_rect_surface(faces['left'],   args.nx, args.ny)
            set_transfinite_rect_surface(faces['right'],  args.nx, args.ny)
            set_transfinite_rect_surface(faces['front'],  args.nx, args.ny)
            set_transfinite_rect_surface(faces['back'],   args.nx, args.ny)
            set_transfinite_rect_surface(faces['bottom'], args.nx, args.ny)
            set_transfinite_rect_surface(faces['top'],    args.nx, args.ny)

        # Generate 3D mesh & write
        gmsh.model.mesh.generate(3)
        gmsh.write(args.outfile)
        print(f"Wrote mesh to {args.outfile}")

        # ---- Verification summary ----
        verify_periodic_nodes(2, [faces['right'], faces['back'], faces['top']])
        verify_periodic_nodes(1, list_curves_on_named_face('right', faces['right']))
        verify_periodic_nodes(0, list_points_on_named_face('right', faces['right']))

    finally:
        gmsh.finalize()

if __name__ == '__main__':
    main()
