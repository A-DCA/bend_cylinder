#!/usr/bin/env python3
"""
Verify periodic boundary conditions in mesh.
"""

import gmsh
import numpy as np
import sys


def verify_periodic_mesh(mesh_file):
    """Verify that periodic boundary conditions were properly applied."""
    
    gmsh.initialize()
    gmsh.open(mesh_file)
    
    print(f"Verifying periodic mesh: {mesh_file}\n")
    
    # Get physical groups
    physical_groups = gmsh.model.getPhysicalGroups()
    inlet_tag = None
    outlet_tag = None
    
    for dim, tag in physical_groups:
        name = gmsh.model.getPhysicalName(dim, tag)
        if name == 'inlet':
            inlet_tag = tag
        elif name == 'outlet':
            outlet_tag = tag
    
    if not inlet_tag or not outlet_tag:
        print("Error: inlet or outlet physical group not found")
        gmsh.finalize()
        return False
    
    # Get inlet and outlet surfaces
    inlet_surfaces = gmsh.model.getEntitiesForPhysicalGroup(2, inlet_tag)
    outlet_surfaces = gmsh.model.getEntitiesForPhysicalGroup(2, outlet_tag)
    
    print(f"Inlet surfaces: {len(inlet_surfaces)}")
    print(f"Outlet surfaces: {len(outlet_surfaces)}")
    
    # Collect all nodes on inlet and outlet
    inlet_coords = {}
    for surf in inlet_surfaces:
        node_tags, coords, _ = gmsh.model.mesh.getNodes(2, surf)
        for i, node_tag in enumerate(node_tags):
            x, y, z = coords[3*i:3*i+3]
            inlet_coords[int(node_tag)] = (x, y, z)
    
    outlet_coords = {}
    for surf in outlet_surfaces:
        node_tags, coords, _ = gmsh.model.mesh.getNodes(2, surf)
        for i, node_tag in enumerate(node_tags):
            x, y, z = coords[3*i:3*i+3]
            outlet_coords[int(node_tag)] = (x, y, z)
    
    print(f"\nNode counts:")
    print(f"  Inlet: {len(inlet_coords)} nodes")
    print(f"  Outlet: {len(outlet_coords)} nodes")
    
    if len(inlet_coords) != len(outlet_coords):
        print("  ✗ MISMATCH: Different number of nodes!")
        gmsh.finalize()
        return False
    else:
        print("  ✓ Same number of nodes")
    
    # Check if each inlet node has a corresponding outlet node
    print("\nChecking node correspondence:")
    
    # Get approximate Z translation
    inlet_z = np.mean([z for x, y, z in inlet_coords.values()])
    outlet_z = np.mean([z for x, y, z in outlet_coords.values()])
    translation_z = outlet_z - inlet_z
    
    print(f"  Inlet Z: {inlet_z:.6f}")
    print(f"  Outlet Z: {outlet_z:.6f}")
    print(f"  Translation Z: {translation_z:.6f}")
    
    # For each inlet node, find matching outlet node
    matched_pairs = 0
    max_distance = 0
    
    inlet_list = list(inlet_coords.values())
    outlet_list = list(outlet_coords.values())
    
    for inlet_pt in inlet_list:
        x_in, y_in, z_in = inlet_pt
        
        # Look for outlet node at same (x,y) position
        best_dist = float('inf')
        for outlet_pt in outlet_list:
            x_out, y_out, z_out = outlet_pt
            
            # Check XY distance (should be near zero)
            dx = x_out - x_in
            dy = y_out - y_in
            dz = z_out - z_in - translation_z
            
            dist = np.sqrt(dx**2 + dy**2 + dz**2)
            if dist < best_dist:
                best_dist = dist
        
        if best_dist < 1e-6:
            matched_pairs += 1
        
        max_distance = max(max_distance, best_dist)
    
    print(f"\nMatching results:")
    print(f"  Matched pairs: {matched_pairs} / {len(inlet_coords)}")
    print(f"  Max XY distance: {max_distance:.2e}")
    
    if matched_pairs == len(inlet_coords) and max_distance < 1e-6:
        print("  ✓ All nodes match - periodic boundary is correct!")
        success = True
    else:
        print("  ✗ Some nodes don't match - periodic boundary may have issues")
        success = False
    
    # Check mesh connectivity at boundaries
    print("\nMesh connectivity check:")
    
    # Get 2D elements on inlet and outlet
    inlet_elem_count = 0
    for surf in inlet_surfaces:
        elem_types, elem_tags, node_tags = gmsh.model.mesh.getElements(2, surf)
        for tags in elem_tags:
            inlet_elem_count += len(tags)
    
    outlet_elem_count = 0
    for surf in outlet_surfaces:
        elem_types, elem_tags, node_tags = gmsh.model.mesh.getElements(2, surf)
        for tags in elem_tags:
            outlet_elem_count += len(tags)
    
    print(f"  Inlet elements: {inlet_elem_count}")
    print(f"  Outlet elements: {outlet_elem_count}")
    
    if inlet_elem_count == outlet_elem_count:
        print("  ✓ Same number of elements")
    else:
        print("  ✗ Different number of elements")
        success = False
    
    gmsh.finalize()
    return success


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python verify_periodic.py <mesh_file.msh>")
        sys.exit(1)
    
    mesh_file = sys.argv[1]
    success = verify_periodic_mesh(mesh_file)
    
    sys.exit(0 if success else 1)
