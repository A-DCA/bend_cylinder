#!/usr/bin/env python3
"""
Sweep annular cross-section (4-arc inner + outer circle) along a curved path.

This creates 5 volumes with proper topology for structured hex meshing:
- 1 inner core volume
- 4 quadrant volumes (annular regions)
"""

import cadquery as cq
import numpy as np
from shapely.geometry import Point
import json
import sys


def circle_intersection(c1_x, c1_y, c2_x, c2_y, r):
    """Find intersection points of two circles with equal radius."""
    d = np.sqrt((c2_x - c1_x)**2 + (c2_y - c1_y)**2)
    if d > 2*r or d < 1e-10:
        return []
    a = d / 2
    h = np.sqrt(r**2 - a**2)
    mx = (c1_x + c2_x) / 2
    my = (c1_y + c2_y) / 2
    dx = (c2_x - c1_x) / d
    dy = (c2_y - c1_y) / d
    p1 = (mx - h * dy, my + h * dx)
    p2 = (mx + h * dy, my - h * dx)
    return [p1, p2]


def is_point_on_boundary(px, py, shape, tolerance=1e-6):
    """Check if a point lies on the boundary of a Shapely polygon."""
    point = Point(px, py)
    return shape.exterior.distance(point) < tolerance


def create_cross_section(radius, offset, outer_radius):
    """Create the annular cross-section geometry."""
    # Calculate 4-circle intersection for inner profile
    circle1 = Point(offset, 0).buffer(radius)
    circle2 = Point(-offset, 0).buffer(radius)
    circle3 = Point(0, offset).buffer(radius)
    circle4 = Point(0, -offset).buffer(radius)
    inner_shape = circle1.intersection(circle2).intersection(circle3).intersection(circle4)
    
    # Find corners
    circle_centers = [(offset, 0), (-offset, 0), (0, offset), (0, -offset)]
    corner_pairs = [(0, 2), (1, 2), (1, 3), (0, 3)]
    
    corners = []
    for i, j in corner_pairs:
        c1 = circle_centers[i]
        c2 = circle_centers[j]
        intersections = circle_intersection(c1[0], c1[1], c2[0], c2[1], radius)
        if intersections:
            for pt in intersections:
                if is_point_on_boundary(pt[0], pt[1], inner_shape, tolerance=radius*0.01):
                    corners.append(pt)
                    break
    
    return corners, circle_centers


def check_inner_profile_fits(corners, outer_radius, radius):
    """
    Check if inner profile fits inside outer circle.
    Returns (is_valid, min_factor, max_distance).
    """
    max_distance = 0.0
    for corner in corners:
        distance = np.sqrt(corner[0]**2 + corner[1]**2)
        max_distance = max(max_distance, distance)
    
    is_valid = max_distance < outer_radius
    min_factor = max_distance / radius
    
    return is_valid, min_factor, max_distance


def create_inner_profile_wire(corners, circle_centers, radius):
    """Create a wire representing the inner 4-arc profile."""
    inner_profile = cq.Workplane("YZ")
    circle_map = [2, 1, 3, 0]
    
    for i in range(4):
        p1 = corners[i]
        p2 = corners[(i + 1) % 4]
        center = circle_centers[circle_map[i]]
        
        angle1 = np.arctan2(p1[1] - center[1], p1[0] - center[0])
        angle2 = np.arctan2(p2[1] - center[1], p2[0] - center[0])
        
        if angle2 - angle1 > np.pi:
            angle2 -= 2*np.pi
        elif angle2 - angle1 < -np.pi:
            angle2 += 2*np.pi
        
        mid_angle = (angle1 + angle2) / 2
        mid_x = center[0] + radius * np.cos(mid_angle)
        mid_y = center[1] + radius * np.sin(mid_angle)
        
        if i == 0:
            inner_profile = inner_profile.moveTo(p1[0], p1[1])
        inner_profile = inner_profile.threePointArc((mid_x, mid_y), (p2[0], p2[1]))
    
    return inner_profile.close()


def create_quadrant_wire(corners, outer_points, circle_centers, radius, outer_radius, i):
    """Create a wire for a single quadrant sector."""
    corner1 = corners[i]
    corner2 = corners[(i + 1) % 4]
    outer1 = outer_points[i]
    outer2 = outer_points[(i + 1) % 4]
    
    sector_sketch = cq.Workplane("YZ")
    sector_sketch = sector_sketch.moveTo(corner1[0], corner1[1])
    sector_sketch = sector_sketch.lineTo(outer1[0], outer1[1])
    
    # Arc along outer circle
    angle1 = np.arctan2(outer1[1], outer1[0])
    angle2 = np.arctan2(outer2[1], outer2[0])
    if angle2 - angle1 > np.pi:
        angle2 -= 2*np.pi
    elif angle2 - angle1 < -np.pi:
        angle2 += 2*np.pi
    mid_angle = (angle1 + angle2) / 2
    mid_x = outer_radius * np.cos(mid_angle)
    mid_y = outer_radius * np.sin(mid_angle)
    
    sector_sketch = sector_sketch.threePointArc((mid_x, mid_y), (outer2[0], outer2[1]))
    sector_sketch = sector_sketch.lineTo(corner2[0], corner2[1])
    
    # Arc along inner profile
    center_idx = [2, 1, 3, 0][i]
    center = circle_centers[center_idx]
    angle1 = np.arctan2(corner2[1] - center[1], corner2[0] - center[0])
    angle2 = np.arctan2(corner1[1] - center[1], corner1[0] - center[0])
    if angle2 - angle1 > np.pi:
        angle2 -= 2*np.pi
    elif angle2 - angle1 < -np.pi:
        angle2 += 2*np.pi
    mid_angle = (angle1 + angle2) / 2
    mid_x = center[0] + radius * np.cos(mid_angle)
    mid_y = center[1] + radius * np.sin(mid_angle)
    
    sector_sketch = sector_sketch.threePointArc((mid_x, mid_y), (corner1[0], corner1[1]))
    return sector_sketch.close()


def create_sinusoidal_path(amplitude, wavelength, n_waves, length, n_points):
    """Create a sinusoidal 3D path as a spline."""
    x = np.linspace(0, length, n_points, endpoint=True)
    z = amplitude * np.sin(2 * np.pi * n_waves * x / length)
    y = np.zeros_like(z)
    
    # Create list of 3D points as tuples
    points = [(float(x[i]), float(y[i]), float(z[i])) for i in range(len(x))]
    
    # Build 3D spline wire using CadQuery
    path = cq.Workplane("XZ").spline(points, forConstruction=True)
    
    return path


def create_straight_path(length):
    """Create a straight vertical path."""
    path = cq.Workplane("XZ").moveTo(0, 0).lineTo(0, length)
    return path


def main():
    if len(sys.argv) < 3:
        print("Usage: python inner_outer_cadquery_curved.py <profile_config.json> <path_config.json>")
        sys.exit(1)
    
    # Load configurations
    profile_config_file = sys.argv[1]
    path_config_file = sys.argv[2]
    
    with open(profile_config_file, 'r') as f:
        profile_config = json.load(f)
    
    with open(path_config_file, 'r') as f:
        path_config = json.load(f)
    
    # Profile parameters
    radius = profile_config['radius']
    offset = profile_config.get('offset')
    if offset is None:
        offset = radius / 1.5
    inner2outer_factor = profile_config.get('inner2outer_factor', 1.5)
    outer_radius = radius * inner2outer_factor
    
    # Get centerline configuration from path config
    centerline = path_config.get('centerline')
    
    if centerline:
        centerline_type = centerline.get('type', 'sinusoidal')
        
        if centerline_type == 'sinusoidal':
            amplitude = centerline['amplitude']
            wavelength = centerline['wavelength']
            n_waves = centerline['n_waves']
            length = centerline['length']
            n_points = centerline.get('n_points', 100)
            
            print(f"Creating curved sweep:")
            print(f"  Centerline: sinusoidal (amplitude={amplitude}, wavelength={wavelength}, n_waves={n_waves})")
            print(f"  Length: {length}")
        elif centerline_type == 'straight':
            length = centerline['length']
            amplitude = 0.0  # No curvature
            
            print(f"Creating straight sweep:")
            print(f"  Length: {length}")
        else:
            print(f"Error: Unknown centerline type: {centerline_type}")
            sys.exit(1)
    else:
        # Straight sweep
        length = path_config.get('length', 0.05)
        amplitude = 0.0
        print(f"Creating straight sweep:")
        print(f"  Length: {length}")
    
    print(f"  Inner radius: {radius}")
    print(f"  Offset: {offset}")
    print(f"  Inner-to-outer factor: {inner2outer_factor}")
    print(f"  Outer radius: {outer_radius}")
    
    # Step 1: Create YZ-plane profile
    print("\nCreating YZ-plane cross-section profile...")
    # Create cross-section geometry
    corners, circle_centers = create_cross_section(radius, offset, outer_radius)
    
    if len(corners) != 4:
        print(f"Error: Expected 4 corners, found {len(corners)}")
        sys.exit(1)
    
    # Check if inner profile fits inside outer circle
    is_valid, min_factor, max_distance = check_inner_profile_fits(corners, outer_radius, radius)
    
    if not is_valid:
        print(f"\n⚠️  Warning: Inner profile does not fit inside outer circle!")
        print(f"   Max distance from origin: {max_distance:.6f}")
        print(f"   Current outer radius: {outer_radius:.6f}")
        print(f"   Current inner2outer_factor: {inner2outer_factor:.3f}")
        print(f"   Minimum required inner2outer_factor: {min_factor:.3f}")
        print(f"\n   Please increase inner2outer_factor to at least {min_factor:.3f}")
        sys.exit(1)
    else:
        print(f"✓ Inner profile validation passed (max distance: {max_distance:.6f} < outer radius: {outer_radius:.6f})")
    
    # Calculate outer points
    outer_points = []
    for corner in corners:
        angle = np.arctan2(corner[1], corner[0])
        outer_x = outer_radius * np.cos(angle)
        outer_y = outer_radius * np.sin(angle)
        outer_points.append((outer_x, outer_y))
    
    # Create profiles
    inner_profile = create_inner_profile_wire(corners, circle_centers, radius)
    
    quadrant_profiles = []
    for i in range(4):
        quadrant_profile = create_quadrant_wire(corners, outer_points, circle_centers, 
                                               radius, outer_radius, i)
        quadrant_profiles.append(quadrant_profile)
    
    # Step 2: Create path
    print("Creating sweep path...")
    if centerline and centerline_type == 'sinusoidal':
        # Curved sweep along spline path
        path = create_sinusoidal_path(amplitude, wavelength, n_waves, length, n_points)
        path_wire = path.val()
        
        print("Sweeping profiles along curved path...")
        inner_core = inner_profile.sweep(path_wire)
        
        quadrants = []
        for i, quadrant_profile in enumerate(quadrant_profiles):
            quadrant = quadrant_profile.sweep(path_wire)
            quadrants.append(quadrant)
        
        output_file = "swept_inner_outer_cylinder_curved.step"
    else:
        # Straight extrusion (for straight centerline or no centerline)
        print("Extruding profiles...")
        inner_core = inner_profile.extrude(length)
        
        quadrants = []
        for quadrant_profile in quadrant_profiles:
            quadrant = quadrant_profile.extrude(length)
            quadrants.append(quadrant)
        
        output_file = "swept_inner_outer_cylinder_straight.step"
    
    # Combine all volumes using add() to preserve separate volumes
    result = inner_core
    for quad in quadrants:
        result = result.add(quad)
    
    cq.exporters.export(result, output_file)
    
    print(f"\nOutput: {output_file}")
    print(f"Corners found: {len(corners)}")


if __name__ == "__main__":
    main()
