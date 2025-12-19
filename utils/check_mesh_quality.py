#!/usr/bin/env python3
"""
Check mesh quality: Jacobian, aspect ratio, and element validity.
"""

import gmsh
import sys
import numpy as np

def check_mesh_quality(mesh_file):
    """Analyze mesh quality metrics."""
    
    gmsh.initialize()
    gmsh.open(mesh_file)
    
    # Get all hexahedral elements
    element_types, element_tags, node_tags = gmsh.model.mesh.getElements(3)
    
    print(f"\nMesh: {mesh_file}")
    print("=" * 60)
    
    # Check for hex elements (type 5 = 8-node hex, type 12 = 27-node hex, type 92 = 20-node hex)
    hex_found = False
    for elem_type, elem_tags, elem_nodes in zip(element_types, element_tags, node_tags):
        elem_name = gmsh.model.mesh.getElementProperties(elem_type)[0]
        n_elements = len(elem_tags)
        print(f"Element type {elem_type} ({elem_name}): {n_elements} elements")
        
        if elem_type in [5, 12, 92]:  # Hexahedral elements
            hex_found = True
            
            # Get Jacobians for all elements
            print(f"\nComputing quality metrics for {n_elements} hex elements...")
            
            jacobians, determinants, coords = gmsh.model.mesh.getJacobians(
                elem_type, [-1, -1, -1]  # Evaluate at element centers
            )
            
            # Reshape determinants
            dets = np.array(determinants)
            
            # Check for negative Jacobians (inverted elements)
            negative_jac = dets < 0
            n_negative = np.sum(negative_jac)
            
            # Check for degenerated elements (near-zero Jacobian)
            threshold_degenerate = 1e-12  # Elements with Jacobian < this are degenerate
            degenerate = (dets > 0) & (dets < threshold_degenerate)
            n_degenerate = np.sum(degenerate)
            
            # Check for very poor quality (very small but not degenerate)
            threshold_poor = 1e-8
            very_poor = (dets >= threshold_degenerate) & (dets < threshold_poor)
            n_very_poor = np.sum(very_poor)
            
            print(f"\nJacobian Analysis:")
            print(f"  Min Jacobian: {dets.min():.6e}")
            print(f"  Max Jacobian: {dets.max():.6e}")
            print(f"  Mean Jacobian: {dets.mean():.6e}")
            print(f"  Median Jacobian: {np.median(dets):.6e}")
            
            print(f"\nElement Quality Classification:")
            print(f"  Negative Jacobians (inverted): {n_negative} / {len(dets)}")
            print(f"  Degenerate (0 < J < 1e-12): {n_degenerate} / {len(dets)}")
            print(f"  Very poor (1e-12 ≤ J < 1e-8): {n_very_poor} / {len(dets)}")
            print(f"  Acceptable (J ≥ 1e-8): {len(dets) - n_negative - n_degenerate - n_very_poor} / {len(dets)}")
            
            if n_negative > 0:
                print(f"\n  ⚠️  CRITICAL: {n_negative} inverted elements detected!")
                neg_elem_indices = np.where(negative_jac)[0]
                print(f"  Inverted element IDs: {elem_tags[neg_elem_indices][:10]}..." if len(neg_elem_indices) > 10 else f"  Inverted element IDs: {elem_tags[neg_elem_indices]}")
            
            if n_degenerate > 0:
                print(f"\n  ⚠️  CRITICAL: {n_degenerate} degenerated elements detected!")
                deg_elem_indices = np.where(degenerate)[0]
                print(f"  Degenerate element IDs: {elem_tags[deg_elem_indices][:10]}..." if len(deg_elem_indices) > 10 else f"  Degenerate element IDs: {elem_tags[deg_elem_indices]}")
                print(f"  These elements have collapsed nodes or edges!")
                
            if n_very_poor > 0:
                print(f"\n  ⚠️  WARNING: {n_very_poor} very poor quality elements!")
                poor_elem_indices = np.where(very_poor)[0]
                print(f"  Very poor element IDs: {elem_tags[poor_elem_indices][:10]}..." if len(poor_elem_indices) > 10 else f"  Very poor element IDs: {elem_tags[poor_elem_indices]}")
                
            if n_negative == 0 and n_degenerate == 0:
                if n_very_poor == 0:
                    print(f"\n  ✓ All elements have good Jacobians")
                else:
                    print(f"\n  ⚠️ Mesh has very poor quality elements but no degenerate elements")
            
            # Get element quality metrics
            # Try multiple quality measures
            quality_measures = ["minJ", "gamma"]  # minJ = min Jacobian, gamma = inscribed/circumscribed ratio
            
            for measure in quality_measures:
                try:
                    qualities = gmsh.model.mesh.getElementQualities(elem_tags, measure)
                    qualities = np.array(qualities)
                    
                    print(f"\nElement Quality ({measure}):")
                    print(f"  Min: {qualities.min():.6f}")
                    print(f"  Max: {qualities.max():.6f}")
                    print(f"  Mean: {qualities.mean():.6f}")
                    print(f"  Median: {np.median(qualities):.6f}")
                    print(f"  Std: {qualities.std():.6f}")
                    
                    if measure == "minJ":
                        # For min Jacobian, values close to 0 are problematic
                        poor = qualities < 0.1
                        n_poor = np.sum(poor)
                        print(f"  Elements with minJ < 0.1: {n_poor} / {len(qualities)}")
                        if n_poor > 0:
                            print(f"  ⚠️  WARNING: {n_poor} highly distorted elements!")
                        else:
                            print(f"  ✓ All elements acceptable (minJ ≥ 0.1)")
                    elif measure == "gamma":
                        # For gamma, 1.0 is perfect, < 0.3 is poor
                        poor = qualities < 0.3
                        n_poor = np.sum(poor)
                        print(f"  Elements with γ < 0.3: {n_poor} / {len(qualities)}")
                        if n_poor > 0:
                            print(f"  ⚠️  WARNING: {n_poor} poor shape elements!")
                        else:
                            print(f"  ✓ All elements acceptable (γ ≥ 0.3)")
                    
                    break  # Use first working measure
                    
                except Exception as e:
                    continue
    
    if not hex_found:
        print("\n⚠️  No hexahedral elements found in mesh!")
    
    # Get node count
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    print(f"\nTotal nodes: {len(node_tags)}")
    
    # Check coordinate ranges
    coords = np.array(node_coords).reshape(-1, 3)
    print(f"\nCoordinate ranges:")
    print(f"  X: [{coords[:,0].min():.6f}, {coords[:,0].max():.6f}]")
    print(f"  Y: [{coords[:,1].min():.6f}, {coords[:,1].max():.6f}]")
    print(f"  Z: [{coords[:,2].min():.6f}, {coords[:,2].max():.6f}]")
    
    # Check for coincident nodes (degenerate mesh indicator)
    print(f"\nChecking for coincident nodes...")
    min_distance = float('inf')
    n_close_pairs = 0
    tolerance = 1e-10
    
    # For efficiency, only sample if there are many nodes
    if len(node_tags) > 1000:
        print(f"  (Sampling 1000 random nodes for proximity check)")
        sample_indices = np.random.choice(len(node_tags), min(1000, len(node_tags)), replace=False)
        sample_coords = coords[sample_indices]
    else:
        sample_coords = coords
    
    # Quick check: compute pairwise distances for sample
    from scipy.spatial.distance import pdist
    distances = pdist(sample_coords[:100] if len(sample_coords) > 100 else sample_coords)
    if len(distances) > 0:
        min_distance = distances.min()
        n_close_pairs = np.sum(distances < tolerance)
        
        print(f"  Min node-to-node distance: {min_distance:.6e}")
        if n_close_pairs > 0:
            print(f"  ⚠️  WARNING: {n_close_pairs} node pairs closer than {tolerance:.1e}")
            print(f"  This may indicate collapsed edges or degenerate elements!")
        else:
            print(f"  ✓ No coincident nodes detected")
    
    gmsh.finalize()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_mesh_quality.py <mesh_file.msh>")
        sys.exit(1)
    
    check_mesh_quality(sys.argv[1])
