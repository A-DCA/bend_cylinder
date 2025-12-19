#!/usr/bin/env python3
import argparse
import gmsh


def face_axes(face_name):
    """Return (u_axis, v_axis) for the face's in-plane directions.
    For axis-aligned unit cube faces:
      - left/right (x=0/1): in-plane axes are y and z
      - front/back (y=0/1): in-plane axes are x and z
      - bottom/top (z=0/1): in-plane axes are x and y
    """
    if face_name in ('left', 'right'):
        return ('y', 'z')
    if face_name in ('front', 'back'):
        return ('x', 'z')
    if face_name in ('bottom', 'top'):
        return ('x', 'y')
    raise ValueError(face_name)


def find_axis_aligned_faces():
    faces = {}
    for dim, tag in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        eps = 1e-6
        if abs(xmin - xmax) < eps:
            if abs(xmin) < eps:
                faces['left'] = tag
            elif abs(xmin - 1.0) < eps or abs(xmax - 1.0) < eps:
                faces['right'] = tag
        elif abs(ymin - ymax) < eps:
            if abs(ymin) < eps:
                faces['front'] = tag
            elif abs(ymin - 1.0) < eps or abs(ymax - 1.0) < eps:
                faces['back'] = tag
        elif abs(zmin - zmax) < eps:
            if abs(zmin) < eps:
                faces['bottom'] = tag
            elif abs(zmin - 1.0) < eps or abs(zmax - 1.0) < eps:
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


def set_transfinite_face(face_name, surf_tag, nx, ny, nz):
    # Decide counts for this face based on its in-plane axes
    u_axis, v_axis = face_axes(face_name)
    axis_to_n = {'x': nx, 'y': ny, 'z': nz}
    edges = get_surface_edges(surf_tag)
    for c in edges:
        ax = edge_axis(c)
        n = axis_to_n[ax]
        gmsh.model.mesh.setTransfiniteCurve(c, n, 'Progression', 1.0)
    gmsh.model.mesh.setTransfiniteSurface(surf_tag)
    gmsh.model.mesh.setRecombine(2, surf_tag)


def set_periodic_surface(target_tag, source_tag, translation_vec):
    tx, ty, tz = translation_vec
    T = [1,0,0,tx,
         0,1,0,ty,
         0,0,1,tz,
         0,0,0,1]
    gmsh.model.mesh.setPeriodic(2, [target_tag], [source_tag], T)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--transfinite', action='store_true')
    ap.add_argument('--nx', type=int, default=11)  # along x-directed edges
    ap.add_argument('--ny', type=int, default=9)   # along y-directed edges
    ap.add_argument('--nz', type=int, default=7)   # along z-directed edges
    ap.add_argument('--outfile', type=str, default='periodic_transfinite_fixed.msh')
    args = ap.parse_args()

    gmsh.initialize()
    try:
        gmsh.model.add('periodic_transfinite_fixed')
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
        gmsh.model.occ.synchronize()

        faces = find_axis_aligned_faces()

        # Mesh sizes to see periodic copy
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.12)
        gmsh.model.mesh.setSize([(0, 1)], 0.02)

        if args.transfinite:
            # Apply axis-aware counts so opposite faces match by axis
            for name in ('left','right','front','back','bottom','top'):
                set_transfinite_face(name, faces[name], args.nx, args.ny, args.nz)

        # Periodic mappings: right<-left (+x), back<-front (+y), top<-bottom (+z)
        set_periodic_surface(faces['right'], faces['left'],  (1,0,0))
        set_periodic_surface(faces['back'],  faces['front'], (0,1,0))
        set_periodic_surface(faces['top'],   faces['bottom'],(0,0,1))

        gmsh.model.mesh.generate(3)
        gmsh.write(args.outfile)
        print(f"Wrote mesh to {args.outfile}")
    finally:
        gmsh.finalize()

if __name__ == '__main__':
    main()
