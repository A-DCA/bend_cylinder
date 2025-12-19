#!/usr/bin/env python3
import argparse
import gmsh


def find_axis_aligned_faces():
    """Return dict with keys 'left','right','front','back','bottom','top' -> surface tags.
    Uses bounding boxes to identify faces of the unit box created at (0,0,0) to (1,1,1).
    """
    faces = {}
    for dim, tag in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        eps = 1e-6
        if abs(xmin - xmax) < eps:  # plane normal along x
            if abs(xmin) < eps:
                faces['left'] = tag
            elif abs(xmin - 1.0) < eps or abs(xmax - 1.0) < eps:
                faces['right'] = tag
        elif abs(ymin - ymax) < eps:  # plane normal along y
            if abs(ymin) < eps:
                faces['front'] = tag  # y=0
            elif abs(ymin - 1.0) < eps or abs(ymax - 1.0) < eps:
                faces['back'] = tag   # y=1
        elif abs(zmin - zmax) < eps:  # plane normal along z
            if abs(zmin) < eps:
                faces['bottom'] = tag  # z=0
            elif abs(zmin - 1.0) < eps or abs(zmax - 1.0) < eps:
                faces['top'] = tag     # z=1
    return faces


def get_surface_edges(surf_tag):
    """Return list of curve tags that bound the given surface, in oriented order."""
    b = gmsh.model.getBoundary([(2, surf_tag)], oriented=True, recursive=False)
    return [t for d, t in b if d == 1]


def edge_axis(curve_tag):
    """Classify edge as 'x','y','z' by looking at its endpoint coordinates."""
    pts = gmsh.model.getBoundary([(1, curve_tag)], oriented=True, recursive=False)
    end_tags = [t for d, t in pts if d == 0]
    xyz = [gmsh.model.getValue(0, p, []) for p in end_tags]
    dx = abs(xyz[0][0] - xyz[1][0])
    dy = abs(xyz[0][1] - xyz[1][1])
    dz = abs(xyz[0][2] - xyz[1][2])
    # Choose axis of largest delta
    if dx >= dy and dx >= dz:
        return 'x'
    elif dy >= dx and dy >= dz:
        return 'y'
    else:
        return 'z'


def set_transfinite_rect_surface(surf_tag, nx, ny):
    """Set transfinite constraints on a rectangular surface: nx along x-edges, ny along y-edges.
    Also set recombine for quads.
    """
    edges = get_surface_edges(surf_tag)
    x_edges = []
    y_edges = []
    for c in edges:
        ax = edge_axis(c)
        if ax == 'x':
            x_edges.append(c)
        elif ax == 'y':
            y_edges.append(c)
        # ignore z-directed edges (shouldn't exist for a planar face)
    # Apply curve constraints
    for c in x_edges:
        gmsh.model.mesh.setTransfiniteCurve(c, nx, 'Progression', 1.0)
    for c in y_edges:
        gmsh.model.mesh.setTransfiniteCurve(c, ny, 'Progression', 1.0)
    # Flag surface
    gmsh.model.mesh.setTransfiniteSurface(surf_tag)
    gmsh.model.mesh.setRecombine(2, surf_tag)


def set_periodic_surface(target_tag, source_tag, translation_vec):
    """Set periodic relation target = copy(source) under translation_vec (length-3)."""
    tx, ty, tz = translation_vec
    T = [1,0,0,tx,
         0,1,0,ty,
         0,0,1,tz,
         0,0,0,1]
    gmsh.model.mesh.setPeriodic(2, [target_tag], [source_tag], T)


def main():
    parser = argparse.ArgumentParser(description='Periodic + Transfinite demo for unit cube')
    parser.add_argument('--transfinite', action='store_true', help='Enable transfinite quad faces with matching counts')
    parser.add_argument('--nx', type=int, default=9, help='#points along x-directed edges of periodic faces')
    parser.add_argument('--ny', type=int, default=7, help='#points along y-directed edges of periodic faces')
    parser.add_argument('--nz', type=int, default=5, help='(unused in 2D faces; for potential extension)')
    parser.add_argument('--outfile', type=str, default='periodic_transfinite_demo.msh', help='Output mesh file')
    args = parser.parse_args()

    gmsh.initialize()
    try:
        gmsh.model.add('periodic_transfinite_demo')
        # Geometry: unit cube
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
        gmsh.model.occ.synchronize()

        # Identify faces by location
        faces = find_axis_aligned_faces()
        # Optionally set local sizes to visualize periodic copy
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.12)
        gmsh.model.mesh.setSize([(0, 1)], 0.02)  # refine one corner point

        
        # Set periodic relations (target <- source) with translations of +1
        set_periodic_surface(faces['right'], faces['left'],  (1,0,0))
        set_periodic_surface(faces['back'],  faces['front'], (0,1,0))
        set_periodic_surface(faces['top'],   faces['bottom'],(0,0,1))

        if args.transfinite:
            # Apply identical transfinite constraints to each periodic pair
            set_transfinite_rect_surface(faces['left'],  args.nx, args.ny)
            set_transfinite_rect_surface(faces['right'], args.nx, args.ny)
            set_transfinite_rect_surface(faces['front'], args.nx, args.ny)
            set_transfinite_rect_surface(faces['back'],  args.nx, args.ny)
            set_transfinite_rect_surface(faces['bottom'], args.nx, args.ny)
            set_transfinite_rect_surface(faces['top'],    args.nx, args.ny)

        # Generate 3D mesh
        gmsh.model.mesh.generate(3)
        gmsh.write(args.outfile)
        print(f"Wrote mesh to {args.outfile}")
    finally:
        gmsh.finalize()


if __name__ == '__main__':
    main()
