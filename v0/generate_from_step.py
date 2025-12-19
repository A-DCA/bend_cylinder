#!/usr/bin/env python3
"""
Generate structured hex mesh from STEP file by adding transfinite constraints.

This script imports a STEP file (e.g., from CadQuery) and applies transfinite
meshing to create a structured hexahedral mesh.
"""

import gmsh
import sys
import json
import argparse


def identify_geometry(step_file):
    """Import STEP file and identify key geometric entities."""
    
    # Initialize GMSH
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    
    # Import the STEP file
    print(f"Importing {step_file}...")
    gmsh.merge(step_file)
    
    # Synchronize to get the model
    gmsh.model.occ.synchronize()
    
    # Get all entities
    volumes = gmsh.model.getEntities(3)
    surfaces = gmsh.model.getEntities(2)
    curves = gmsh.model.getEntities(1)
    points = gmsh.model.getEntities(0)
    
    print(f"\nGeometry entities:")
    print(f"  Volumes: {len(volumes)}")
    print(f"  Surfaces: {len(surfaces)}")
    print(f"  Curves: {len(curves)}")
    print(f"  Points: {len(points)}")
    
    # Get bounding boxes to identify surfaces
    surface_info = []
    for dim, tag in surfaces:
        bbox = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox
        center = ((xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2)
        
        # Classify surface by position and geometry
        surface_type = "unknown"
        length_z = zmax - zmin
        
        if abs(length_z) < 1e-6:  # Flat in Z (inlet/outlet)
            if abs(zmin) < 1e-6:
                surface_type = "inlet"
            else:
                surface_type = "outlet"
        else:
            # Check if it's a curved outer wall or flat radial plane
            # Outer walls are cylindrical surfaces far from origin
            r_center = (center[0]**2 + center[1]**2)**0.5
            r_min = (xmin**2 + ymin**2)**0.5
            r_max = (xmax**2 + ymax**2)**0.5
            
            if abs(r_max - r_min) < 1e-6 and r_center > 0.006:  # Cylindrical outer wall
                surface_type = "wall"
            elif abs(r_max - r_min) > 1e-6:  # Radial variation
                surface_type = "radial_or_arc"
        
        surface_info.append({
            'tag': tag,
            'bbox': bbox,
            'center': center,
            'type': surface_type
        })
    
    return volumes, surfaces, curves, points, surface_info


def apply_transfinite_meshing(n_vertical, n_circumferential, n_radial, surface_info):
    """Apply transfinite constraints to curves and suitable surfaces."""
    
    curves = gmsh.model.getEntities(1)
    surfaces = gmsh.model.getEntities(2)
    volumes = gmsh.model.getEntities(3)
    
    print(f"\nApplying transfinite constraints:")
    print(f"  n_vertical = {n_vertical}")
    print(f"  n_circumferential = {n_circumferential}")
    print(f"  n_radial = {n_radial}")
    
    # Classify curves by orientation and apply transfinite
    vertical_curves = []
    circumferential_curves = []
    radial_curves = []
    
    # Debug: collect curve types
    curve_types_found = set()
    
    for dim, tag in curves:
        bbox = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox
        
        length_x = xmax - xmin
        length_y = ymax - ymin
        length_z = zmax - zmin
        max_xy = max(length_x, length_y)
        
        # Get curve type to distinguish straight lines from arcs
        curve_type = gmsh.model.getType(dim, tag)
        curve_types_found.add(curve_type)
        
        # Calculate chord length (straight distance between endpoints)
        # For swept geometries, use arc/chord ratio to detect straightness
        chord_length = ((xmax - xmin)**2 + (ymax - ymin)**2 + (zmax - zmin)**2)**0.5
        
        # Sample points along curve to estimate arc length
        try:
            # Get parametric bounds
            bounds = gmsh.model.getParametrizationBounds(dim, tag)
            t_min, t_max = bounds[0][0], bounds[1][0]
            n_samples = 20
            arc_length = 0.0
            prev_pt = gmsh.model.getValue(dim, tag, [t_min])
            for i in range(1, n_samples + 1):
                t = t_min + i * (t_max - t_min) / n_samples
                pt = gmsh.model.getValue(dim, tag, [t])
                arc_length += sum((pt[j] - prev_pt[j])**2 for j in range(3))**0.5
                prev_pt = pt
            curvature_ratio = arc_length / chord_length if chord_length > 1e-10 else 1.0
        except:
            curvature_ratio = 1.0
        
        # Classify by geometry and orientation
        if length_z < 1e-6 and max_xy > 1e-6:
            # Flat in Z - in XY plane (inlet/outlet cross-sections)
            if curvature_ratio > 1.02:
                # Curved line - circumferential arc
                gmsh.model.mesh.setTransfiniteCurve(tag, n_circumferential)
                circumferential_curves.append(tag)
            else:
                # Nearly straight - radial line at inlet/outlet
                gmsh.model.mesh.setTransfiniteCurve(tag, n_radial)
                radial_curves.append(tag)
        elif length_z > max_xy * 5:
            # Strongly Z-dominant - vertical curves along sweep path
            gmsh.model.mesh.setTransfiniteCurve(tag, n_vertical)
            vertical_curves.append(tag)
        else:
            # Mixed Z and XY - distinguish by curvature
            if curvature_ratio > 1.02:
                # Curved - swept circumferential arc
                gmsh.model.mesh.setTransfiniteCurve(tag, n_circumferential)
                circumferential_curves.append(tag)
            else:
                # Nearly straight - swept radial line
                gmsh.model.mesh.setTransfiniteCurve(tag, n_radial)
                radial_curves.append(tag)
    
    print(f"  Curve types found: {curve_types_found}")
    print(f"  Classified {len(vertical_curves)} vertical curves")
    print(f"  Classified {len(circumferential_curves)} circumferential curves")
    print(f"  Classified {len(radial_curves)} radial curves")
    
    # Apply transfinite to surfaces with 4 boundary curves
    transfinite_surfaces = []
    for s in surface_info:
        tag = s['tag']
        stype = s['type']
        
        boundary = gmsh.model.getBoundary([(2, tag)], oriented=False, recursive=False)
        n_boundary_curves = len(boundary)
        
        # Try transfinite on 4-sided surfaces (ruled surfaces)
        if n_boundary_curves == 4:
            try:
                gmsh.model.mesh.setTransfiniteSurface(tag)
                gmsh.model.mesh.setRecombine(2, tag)
                transfinite_surfaces.append(tag)
            except Exception as e:
                # If transfinite fails, still recombine
                gmsh.model.mesh.setRecombine(2, tag)
        else:
            # Recombine for quad mesh where possible
            gmsh.model.mesh.setRecombine(2, tag)
    
    print(f"  Applied transfinite to {len(transfinite_surfaces)} surfaces")
    
    # Try to apply transfinite to volumes with exactly 6 faces
    transfinite_volumes = []
    for dim, tag in volumes:
        boundary = gmsh.model.getBoundary([(3, tag)], oriented=False, recursive=False)
        n_faces = len(boundary)
        
        if n_faces == 6:
            try:
                gmsh.model.mesh.setTransfiniteVolume(tag)
                transfinite_volumes.append(tag)
                print(f"  Volume {tag}: transfinite (6 faces)")
            except Exception as e:
                print(f"  Volume {tag}: skipped - {e}")
        else:
            print(f"  Volume {tag}: {n_faces} faces (needs 6 for transfinite)")
    
    if transfinite_volumes:
        print(f"\n  Applied transfinite to {len(transfinite_volumes)} volumes (will generate hex mesh)")
    else:
        print(f"\n  Note: No volumes suitable for transfinite (will use tetrahedral mesh)")


def check_periodic_compatibility(inlet_surfaces, outlet_surfaces):
    """Check if inlet and outlet surfaces can form periodic boundaries."""
    
    if not inlet_surfaces or not outlet_surfaces:
        print("  Cannot check periodicity: missing inlet or outlet")
        return False, None
    
    print("\n  Checking periodic boundary compatibility...")
    
    # For each inlet surface, find matching outlet surface
    periodic_pairs = []
    
    for inlet_tag in inlet_surfaces:
        # Get inlet surface center and area
        inlet_bbox = gmsh.model.getBoundingBox(2, inlet_tag)
        inlet_center_x = (inlet_bbox[0] + inlet_bbox[3]) / 2
        inlet_center_y = (inlet_bbox[1] + inlet_bbox[4]) / 2
        inlet_z = inlet_bbox[2]  # Should be near 0
        
        # Find outlet surface with matching XY position
        for outlet_tag in outlet_surfaces:
            outlet_bbox = gmsh.model.getBoundingBox(2, outlet_tag)
            outlet_center_x = (outlet_bbox[0] + outlet_bbox[3]) / 2
            outlet_center_y = (outlet_bbox[1] + outlet_bbox[4]) / 2
            outlet_z = outlet_bbox[5]  # Should be at max Z
            
            # Check if XY centers match (same annular cross-section)
            dx = abs(outlet_center_x - inlet_center_x)
            dy = abs(outlet_center_y - inlet_center_y)
            
            if dx < 1e-6 and dy < 1e-6:
                # Found matching pair
                length = outlet_z - inlet_z
                periodic_pairs.append({
                    'source': inlet_tag,
                    'target': outlet_tag,
                    'translation': [0, 0, length]
                })
                print(f"    Surface pair ({inlet_tag}, {outlet_tag}): compatible for periodicity")
                print(f"      Translation: [0, 0, {length:.6f}]")
    
    if periodic_pairs:
        print(f"  Found {len(periodic_pairs)} periodic surface pair(s)")
        return True, periodic_pairs
    else:
        print("  No matching periodic surface pairs found")
        return False, None


def apply_periodic_constraints(periodic_pairs):
    """Apply periodic boundary constraints to matched surface pairs."""
    
    if not periodic_pairs:
        return
    
    print("\n  Applying periodic constraints...")
    
    # Increase the tolerance for matching points significantly
    # The geometry has some numerical precision issues from STEP import
    gmsh.option.setNumber("Geometry.Tolerance", 1e-5)
    gmsh.option.setNumber("Geometry.ToleranceBoolean", 1e-5)
    gmsh.option.setNumber("Mesh.ToleranceInitialDelaunay", 1e-5)
    
    for pair in periodic_pairs:
        source_tag = pair['source']
        target_tag = pair['target']
        translation = pair['translation']
        
        try:
            # Define affine transformation: translation only
            # GMSH expects a 4x4 transformation matrix (16 values) in row-major format
            affine_transform = [
                1, 0, 0, translation[0],  # Row 1
                0, 1, 0, translation[1],  # Row 2
                0, 0, 1, translation[2],  # Row 3
                0, 0, 0, 1                # Row 4
            ]
            
            # Set periodic constraint (this also handles mesh periodicity automatically)
            gmsh.model.mesh.setPeriodic(2, [target_tag], [source_tag], affine_transform)
            print(f"    ✓ Applied periodic constraint: surface {target_tag} → {source_tag}")
            
        except Exception as e:
            print(f"    ✗ Failed to apply periodic constraint ({source_tag}, {target_tag}): {e}")


def define_physical_groups(surface_info, enable_periodic=False):
    """Define physical groups for boundaries."""
    
    inlet_surfaces = [s['tag'] for s in surface_info if s['type'] == 'inlet']
    outlet_surfaces = [s['tag'] for s in surface_info if s['type'] == 'outlet']
    wall_surfaces = [s['tag'] for s in surface_info if s['type'] == 'wall']
    
    # If no walls detected by type, find cylindrical outer surfaces
    if not wall_surfaces:
        for s in surface_info:
            if s['type'] == 'radial_or_arc':
                # Check if it's far from origin (outer wall)
                bbox = s['bbox']
                xmin, ymin, zmin, xmax, ymax, zmax = bbox
                r_max = max((xmin**2 + ymin**2)**0.5, (xmax**2 + ymax**2)**0.5)
                if r_max > 0.006:  # Outer radius region
                    wall_surfaces.append(s['tag'])
    
    volumes = gmsh.model.getEntities(3)
    volume_tags = [tag for dim, tag in volumes]
    
    # Check and apply periodic boundaries if requested
    periodic_pairs = None
    if enable_periodic:
        is_periodic, periodic_pairs = check_periodic_compatibility(inlet_surfaces, outlet_surfaces)
        if is_periodic:
            # Note: periodic constraints are applied BEFORE mesh generation
            # We'll return the pairs to apply them at the right time
            pass
    
    if inlet_surfaces:
        gmsh.model.addPhysicalGroup(2, inlet_surfaces, name="inlet")
        print(f"  inlet: surfaces {inlet_surfaces}")
    
    if outlet_surfaces:
        gmsh.model.addPhysicalGroup(2, outlet_surfaces, name="outlet")
        print(f"  outlet: surfaces {outlet_surfaces}")
    
    if wall_surfaces:
        gmsh.model.addPhysicalGroup(2, wall_surfaces, name="walls")
        print(f"  walls: surfaces {wall_surfaces}")
    
    if volume_tags:
        gmsh.model.addPhysicalGroup(3, volume_tags, name="fluid")
        print(f"  fluid: volumes {volume_tags}")
    
    return periodic_pairs


def generate_mesh(step_file, output_file, n_vertical=21, n_circumferential=20, n_radial=10, periodic=False):
    """Generate structured hex mesh from STEP file."""
    
    # Identify geometry
    volumes, surfaces, curves, points, surface_info = identify_geometry(step_file)
    
    # Apply transfinite meshing
    apply_transfinite_meshing(n_vertical, n_circumferential, n_radial, surface_info)
    
    # Define physical groups and check periodicity
    print("\nDefining physical groups:")
    periodic_pairs = define_physical_groups(surface_info, enable_periodic=periodic)
    
    # Apply periodic constraints BEFORE mesh generation
    if periodic and periodic_pairs:
        apply_periodic_constraints(periodic_pairs)
    
    # Generate mesh
    print("\nGenerating mesh...")
    gmsh.model.mesh.generate(3)
    
    # Set mesh order to 2 for curved elements
    print("Creating high-order curved mesh...")
    gmsh.model.mesh.setOrder(2)
    
    # Get mesh statistics
    nodes = gmsh.model.mesh.getNodes()
    elements = gmsh.model.mesh.getElements()
    n_nodes = len(nodes[0])
    n_elements = sum(len(e) for e in elements[1])
    
    print(f"\nMesh statistics:")
    print(f"  Nodes: {n_nodes}")
    print(f"  Elements: {n_elements}")
    
    # Save mesh
    gmsh.write(output_file)
    print(f"\nMesh written to: {output_file}")
    
    # Cleanup
    gmsh.finalize()


def main():
    parser = argparse.ArgumentParser(
        description="Generate structured hex mesh from STEP file"
    )
    parser.add_argument("step_file", help="Input STEP file")
    parser.add_argument("output_file", help="Output mesh file (.msh)")
    parser.add_argument("--n-vertical", type=int, default=21,
                       help="Number of divisions in vertical direction")
    parser.add_argument("--n-circumferential", type=int, default=20,
                       help="Number of divisions in circumferential direction")
    parser.add_argument("--n-radial", type=int, default=10,
                       help="Number of divisions in radial direction")
    parser.add_argument("--periodic", action="store_true",
                       help="Enable periodic boundary conditions between inlet and outlet")
    
    args = parser.parse_args()
    
    generate_mesh(
        args.step_file,
        args.output_file,
        args.n_vertical,
        args.n_circumferential,
        args.n_radial,
        args.periodic
    )


if __name__ == "__main__":
    main()
