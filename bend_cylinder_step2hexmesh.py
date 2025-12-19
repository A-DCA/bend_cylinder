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

#types of volumes
#inner_core - 0
#quadrant0 -1
#quadrant1 -2
#quadrant2 -3
#quadrant -4

def classify_volume_types(center_points_inlet_plane, epsilon=1.0e-8):
    x,y, z = center_points_inlet_plane
    if(abs(y) < epsilon and abs(z) < epsilon):
        return 0 #"inner_core"
    
    if(abs(y) < epsilon and z < 0 and  abs(z) > epsilon):
        return 1 #"quadrant0"
    
    if(y > epsilon and abs(z) < epsilon):
        return 2 #"quadrant1"
    
    if(abs(y) < epsilon and z > epsilon):
            return 3 #"quadrant2"

    if(y <0 and abs(y) > epsilon and abs(z) < epsilon):
        return 4 #"quadrant3"

#types of surface
#inlet_plane
#outlet_plane
#radial_plane
#curved_outer_wall
#curved_inner_wall

#types of edges in inlet_plane
#radial_edge
#outer_circumferential_edge   
#inner_circumferential_edge
#along_edge

def edg_length(edge_tag, n_samples=20):
    bounds = gmsh.model.getParametrizationBounds(1, edge_tag)
    t_min, t_max = bounds[0][0], bounds[1][0]
    arc_length = 0.0
    prev_pt = gmsh.model.getValue(1, edge_tag, [t_min])
    for i in range(1, n_samples + 1):
        t = t_min + i * (t_max - t_min) / n_samples
        pt = gmsh.model.getValue(1, edge_tag, [t])
        arc_length += sum((pt[j] - prev_pt[j])**2 for j in range(3))**0.5
        prev_pt = pt
    return arc_length

def isEqual(edge_tags, i, epsilon=1.0e-8):
    edge_tag1 = edge_tags[i]
    edge_tag2 = edge_tags[(i+2) % 4]

    if abs(edg_length(edge_tag1, n_samples=20) - edg_length(edge_tag2, n_samples=20)) > epsilon:
        return False    
    else:
        return True


def label_edges(edge_tags, volume_label, epsilon=1.0e-8):
    if(len(edge_tags) !=4):
        raise ValueError(f"Edge tags length is not 4.")

    edge_types = {}
    for i in range(4):
        edge_tag0 = edge_tags[i]
        edge_tag2 = edge_tags[(i+2) % 4]

        edge_tag1 = edge_tags[(i+1) % 4]
        edge_tag3 = edge_tags[(i+3) % 4]
        if volume_label == 0:  
            if isEqual(edge_tags, i, epsilon):
                edge_types[edge_tag0] = "inner_circumferential_edge"
                edge_types[edge_tag2] = "inner_circumferential_edge"
                edge_types[edge_tag1] = "inner_circumferential_edge"
                edge_types[edge_tag2] = "inner_circumferential_edge"
        else:
            if isEqual(edge_tags, i, epsilon):
                edge_types[edge_tag0] = "radial_edge"    
                edge_types[edge_tag2] = "radial_edge"  
                if  edg_length(edge_tag1, n_samples=20) < edg_length(edge_tag3, n_samples=20):
                    edge_types[edge_tag1] = "inner_circumferential_edge"    
                    edge_types[edge_tag3] = "outer_circumferential_edge" 
                else:
                    edge_types[edge_tag3] = "inner_circumferential_edge"    
                    edge_types[edge_tag1] = "outer_circumferential_edge"
            elif isEqual(edge_tags, (i+1) %4, epsilon):
                edge_types[edge_tag1] = "radial_edge"    
                edge_types[edge_tag3] = "radial_edge"    
                if  edg_length(edge_tag0, n_samples=20) < edg_length(edge_tag2, n_samples=20):
                    edge_types[edge_tag0] = "inner_circumferential_edge"    
                    edge_types[edge_tag2] = "outer_circumferential_edge" 
                else:
                    edge_types[edge_tag2] = "inner_circumferential_edge"    
                    edge_types[edge_tag0] = "outer_circumferential_edge"

    return edge_types

def is_rotation_tuples(lst1, lst2, tol=1e-8):
    if len(lst1) != len(lst2):
        return False
    n = len(lst1)
    for k in range(n):
        if all(all(abs(a - b) < tol for a, b in zip(lst1[(i + k) % n], lst2[i])) for i in range(n)):
            return True
    return False

def match_inlet_outlet(inlet_plane_edge_types, outlet_plane_edge_types, epsilon=1.0e-8):
    matched = [False] * len(inlet_plane_edge_types)
    offset = [-1] * len(inlet_plane_edge_types)
    for n, (edge_tag1, edge_type) in enumerate(inlet_plane_edge_types.items()):
        x1 = 0.0
        edge_tag1_len = edg_length(edge_tag1, n_samples=20)
        point_tags1 = gmsh.model.getBoundary([(1, edge_tag1)],             
                                        combined=False, oriented=False, recursive=False)
        coordinates1 = []
        for dim, pt_tag in point_tags1:
            x,y,z = gmsh.model.getValue(0, pt_tag, [])
            x1 += x
            coordinates1.append((y, z))
        x1 /= len(point_tags1)

        print(f"Matching inlet edge1 {edge_tag1} of type {edge_type} and length {edge_tag1_len:.6f}")
        print(f"  Coordinates1: {coordinates1}")
        for edge_tag2, edge_type2 in outlet_plane_edge_types.items():
            if edge_type == edge_type2:
                edge_tag2_len = edg_length(edge_tag2, n_samples=20)
                point_tags2 = gmsh.model.getBoundary([(1, edge_tag2)], 
                                        combined=False, oriented=False, recursive=False)      
                coordinates2 = []
                x2 = 0.0

                for dim, pt_tag in point_tags2:
                    x,y,z = gmsh.model.getValue(0, pt_tag, [])
                    coordinates2.append((y, z))
                    x2 += x
                x2 /= len(point_tags2)

                print(f"Matching inlet edge2 {edge_tag2} of type {edge_type2} and length {edge_tag2_len:.6f}")
                print(f"  Coordinates2: {coordinates2}")
                if abs(edge_tag1_len - edge_tag2_len) < epsilon:     
                    matched[n] = is_rotation_tuples(coordinates1, coordinates2, tol=epsilon)
                    if matched[n]:
                        offset[n] = x2 - x1
                        break
    return matched[0] and matched[1] and matched[2] and matched[3], (offset[0]+offset[1]+offset[2]+offset[3])/4 

def label_entities(epsilon=1.0e-8):
    print(f"{len(gmsh.model.getEntities(3))} volumes to be processed.")
    print(f"{len(gmsh.model.getEntities(2))} surfaces to be processed.")
    print(f"{len(gmsh.model.getEntities(1))} edges to be processed.")
    print(f"{len(gmsh.model.getEntities(0))} vertices to be processed.")

    volume_info = {}
    surface_info = {}
    edge_info = {}

    epsilon = 1e-6
    for dim, vol_tag in gmsh.model.getEntities(3):  # volumes
        volume_info[vol_tag] = (-1, False)  # (type, isPeriodic)
        edge_info[vol_tag] = {}
        surface_info[vol_tag] = {}

        surf_tags = gmsh.model.getBoundary([(dim, vol_tag)], 
                                    combined=False, oriented=False, recursive=False)
        print(f"Volume {vol_tag} has {len(surf_tags)} surfaces: {surf_tags}.")
        for dim, surf_tag  in surf_tags: 
            surface_info[vol_tag][surf_tag] = "unknown"

            edge_info[vol_tag][surf_tag] = {}

            edge_tags = gmsh.model.getBoundary([(dim, surf_tag)], 
                                        combined=False, oriented=False, recursive=False)
            print(f"Surface {surf_tag} has {len(edge_tags)} edges: {edge_tags}.")
            for dim, edge_tag in edge_tags:
                edge_info[vol_tag][surf_tag][edge_tag] = "along_edge"
                vertex_tags = gmsh.model.getBoundary([(dim, edge_tag)],
                                        combined=False, oriented=False, recursive=False) 
                print(f"Edge {edge_tag} has {len(vertex_tags)} vertices: {vertex_tags}.")


    for dim, vol_tag in gmsh.model.getEntities(3):  # volumes
        surf_tags = gmsh.model.getBoundary([(dim, vol_tag)], 
                                    combined=False, oriented=False, recursive=False)
        inlet_plane_edges = (0,[])
        outlet_plane_edges = (0,[])
        for dim, surf_tag in surf_tags:
            surf_type = gmsh.model.getType(dim, surf_tag)
            if(surf_type == "Plane"):
                u, v = 0.5, 0.5
                nx, ny, nz = gmsh.model.getNormal(surf_tag, [u, v])
                if abs(nx+1) <epsilon and abs(ny) <epsilon and abs(nz) <epsilon:
                    surf_type = "inlet_plane"                    
                    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, surf_tag)
                    center = [(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2]
                    volume_info[vol_tag] = (classify_volume_types(center), False)
                    edge_tags = gmsh.model.getBoundary([(dim, surf_tag)], 
                                        combined=False, oriented=False, recursive=False)
                    inlet_plane_edges = (surf_tag, edge_tags)
                elif abs(nx-1) <epsilon and abs(ny) <epsilon and abs(nz) <epsilon:
                    surf_type = "outlet_plane"        
                    edge_tags = gmsh.model.getBoundary([(dim, surf_tag)], 
                                        combined=False, oriented=False, recursive=False)
                    outlet_plane_edges = (surf_tag, edge_tags)          
            surface_info[vol_tag][surf_tag] = surf_type
        
        print(f"Volume {vol_tag} info: {volume_info[vol_tag]}: surfaces {surf_tags} classified as type: \033[91m{surface_info[vol_tag]}\033[0m.")
            
        #print(f"Edges in surface {surf_tag}: {edge_info[vol_tag][surf_tag]}")
        
        print("inlet_plane_edges:", inlet_plane_edges)
        print("outlet_plane_edges:", outlet_plane_edges)

        surf_tag, inlet_edges = inlet_plane_edges
        if(len(inlet_edges) !=4):
            raise ValueError(f"Inlet plane of volume {vol_tag} does not have 4 edges.")
        #same length inlet_plane_edges[i] and inlet_plane_edges[(i+2) % 4]

        inlet_plane_edge_types = {}
        edge_types = label_edges([edge_tag for dim, edge_tag in inlet_edges], volume_info[vol_tag][0], epsilon)    
        for edge_tag, edge_type in edge_types.items():
            edge_info[vol_tag][surf_tag][edge_tag] = edge_type
            inlet_plane_edge_types[edge_tag] = edge_type
            up, _ = gmsh.model.getAdjacencies(1, edge_tag)
            for surf_tag2 in up:
                if surf_tag2 != surf_tag and surf_tag2 in surface_info[vol_tag]:
                    edge_info[vol_tag][surf_tag2][edge_tag] = edge_type
    

            surf_tag, outlet_edges = outlet_plane_edges
        if(len(outlet_edges) !=4):
            raise ValueError(f"Outlet plane of volume {vol_tag} does not have 4 edges.")
        #same length inlet_plane_edges[i] and inlet_plane_edges[(i+2) % 4]

        outlet_plane_edge_types = {}
        edge_types = label_edges([edge_tag for dim, edge_tag in outlet_edges], volume_info[vol_tag][0], epsilon)    
        for edge_tag, edge_type in edge_types.items():
            edge_info[vol_tag][surf_tag][edge_tag] = edge_type
            outlet_plane_edge_types[edge_tag] = edge_type
            up, _ = gmsh.model.getAdjacencies(1, edge_tag)
            for surf_tag2 in up:
                if surf_tag2 != surf_tag and surf_tag2 in surface_info[vol_tag]:
                    edge_info[vol_tag][surf_tag2][edge_tag] = edge_type


        for dim, surf_tag in surf_tags:
            if(surface_info[vol_tag][surf_tag] != "outlet_plane" and 
               surface_info[vol_tag][surf_tag] != "inlet_plane"):
                edge_tags = gmsh.model.getBoundary([(dim, surf_tag)], 
                                        combined=False, oriented=False, recursive=False)
                common_edges = set(inlet_edges).intersection(set(edge_tags))
                for edge_dim, edge_tag in common_edges:
                    if(volume_info[vol_tag][0] == 0):
                        if(inlet_plane_edge_types[edge_tag] == "inner_circumferential_edge"):
                            surface_info[vol_tag][surf_tag] = "inner_circumferential_wall"
                    else:
                        if(inlet_plane_edge_types[edge_tag] == "inner_circumferential_edge"):
                            surface_info[vol_tag][surf_tag] = "inner_circumferential_wall"
                        if(inlet_plane_edge_types[edge_tag] ==  "outer_circumferential_edge"):
                            surface_info[vol_tag][surf_tag] = "outer_circumferential_wall"
                        if(inlet_plane_edge_types[edge_tag] == "radial_edge"):
                            surface_info[vol_tag][surf_tag] = "curved_outer_wall"
                        if(inlet_plane_edge_types[edge_tag] == "radial_edge"):
                            surface_info[vol_tag][surf_tag] = "radial_plane"

        periodic = match_inlet_outlet(inlet_plane_edge_types, outlet_plane_edge_types)
        
        volume_info[vol_tag] = (volume_info[vol_tag][0], periodic)


    for dim, vol_tag in  gmsh.model.getEntities(3):
        print(f"\nVolume {vol_tag} classified as type {volume_info[vol_tag][0]}")
        surface_tags = gmsh.model.getBoundary([(dim, vol_tag)], 
                                        combined=False, oriented=False, recursive=False)
        for dim, surf_tag in surface_tags:    
            print(f"surface_info[{vol_tag}][{surf_tag}] ={surface_info[vol_tag][surf_tag]}")        
            edge_tags = gmsh.model.getBoundary([(dim, surf_tag)], 
                                        combined=False, oriented=False, recursive=False)
            for dim, edge_tag in edge_tags:    
                print(f"  edge_info[{vol_tag}][{surf_tag}][{edge_tag}] = {edge_info[vol_tag][surf_tag][edge_tag]}")
    
    return volume_info, surface_info, edge_info

def label_entities0(epsilon=1.0e-8):
    
    #lable volumes according to inlet plane mass center locations, five points on the same planes 
    #  
    volumes = gmsh.model.getEntities(3)
    volume_info = {}
    surface_info = {}
    edge_info = {}

    for dim, vol_tag in volumes:
        surface_info[vol_tag] = {}
        surface_tags = gmsh.model.getBoundary([(dim, vol_tag)], 
                                        combined=False, oriented=False, recursive=False)
        if(len(surface_tags)!=6):
            raise ValueError(f"Volume {vol_tag} does not have 6 boundary surfaces.")
        
        edge_info[vol_tag] = {}
        for dim, surf_tag in surface_tags:
            edge_info[vol_tag][surf_tag] = {}

            edge_tags = gmsh.model.getBoundary([(dim, surf_tag)], 
                                        combined=False, oriented=False, recursive=False)
            for dim, edge_tag in edge_tags:
                edge_info[vol_tag][surf_tag][edge_tag] = "along_edge"
            
        inlet_plane_edges = (0,[])
        outlet_plane_edges = (0,[])

        for dim, surf_tag in surface_tags:
            surf_type = gmsh.model.getType(dim, surf_tag)
            if(surf_type == "Plane"):
                u, v = 0.5, 0.5
                nx, ny, nz = gmsh.model.getNormal(surf_tag, [u, v])
                if abs(nx+1) <epsilon and abs(ny) <epsilon and abs(nz) <epsilon:
                    surf_type = "inlet_plane"                    
                    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, surf_tag)
                    center = [(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2]
                    volume_info[vol_tag] = (classify_volume_types(center), False)
                    edge_tags = gmsh.model.getBoundary([(dim, surf_tag)], 
                                        combined=False, oriented=False, recursive=False)
                    inlet_plane_edges = (surf_tag, edge_tags)
                elif abs(nx-1) <epsilon and abs(ny) <epsilon and abs(nz) <epsilon:
                    surf_type = "outlet_plane"        
                    edge_tags = gmsh.model.getBoundary([(dim, surf_tag)], 
                                        combined=False, oriented=False, recursive=False)
                    outlet_plane_edges = (surf_tag, edge_tags)          
            surface_info[vol_tag][surf_tag] = surf_type


        surf_tag, inlet_edges = inlet_plane_edges
        if(len(inlet_edges) !=4):
            raise ValueError(f"Inlet plane of volume {vol_tag} does not have 4 edges.")
        #same length inlet_plane_edges[i] and inlet_plane_edges[(i+2) % 4]

        inlet_plane_edge_types = {}
        edge_types = label_edges([edge_tag for dim, edge_tag in inlet_edges], volume_info[vol_tag][0], epsilon)    
        for edge_tag, edge_type in edge_types.items():
            edge_info[vol_tag][surf_tag][edge_tag] = edge_type
            inlet_plane_edge_types[edge_tag] = edge_type
            up, _ = gmsh.model.getAdjacencies(1, edge_tag)
            for surf_tag2 in up:
                if surf_tag2 != surf_tag:
                    edge_info[vol_tag][surf_tag2][edge_tag] = edge_type

        surf_tag, outlet_edges = outlet_plane_edges
        if(len(outlet_edges) !=4):
            raise ValueError(f"Outlet plane of volume {vol_tag} does not have 4 edges.")
        #same length inlet_plane_edges[i] and inlet_plane_edges[(i+2) % 4]

        outlet_plane_edge_types = {}
        edge_types = label_edges([edge_tag for dim, edge_tag in outlet_edges], volume_info[vol_tag][0], epsilon)    
        for edge_tag, edge_type in edge_types.items():
            edge_info[vol_tag][surf_tag][edge_tag] = edge_type
            outlet_plane_edge_types[edge_tag] = edge_type
            up, _ = gmsh.model.getAdjacencies(1, edge_tag)
            for surf_tag2 in up:
                if surf_tag2 != surf_tag:
                    edge_info[vol_tag][surf_tag2][edge_tag] = edge_type


        for dim, surf_tag in surface_tags:
            if(surface_info[vol_tag][surf_tag] != "outlet_plane" and 
               surface_info[vol_tag][surf_tag] != "inlet_plane"):
                edge_tags = gmsh.model.getBoundary([(dim, surf_tag)], 
                                        combined=False, oriented=False, recursive=False)
                common_edges = set(inlet_edges).intersection(set(edge_tags))
                for edge_dim, edge_tag in common_edges:
                    if(volume_info[vol_tag][0] == 0):
                        if(inlet_plane_edge_types[edge_tag] == "inner_circumferential_edge"):
                            surface_info[vol_tag][surf_tag] = "inner_circumferential_wall"
                    else:
                        if(inlet_plane_edge_types[edge_tag] == "inner_circumferential_edge"):
                            surface_info[vol_tag][surf_tag] = "inner_circumferential_wall"
                        if(inlet_plane_edge_types[edge_tag] ==  "outer_circumferential_edge"):
                            surface_info[vol_tag][surf_tag] = "outer_circumferential_wall"
                        if(inlet_plane_edge_types[edge_tag] == "radial_edge"):
                            surface_info[vol_tag][surf_tag] = "curved_outer_wall"
                        if(inlet_plane_edge_types[edge_tag] == "radial_edge"):
                            surface_info[vol_tag][surf_tag] = "radial_plane"

        periodic = match_inlet_outlet(inlet_plane_edge_types, outlet_plane_edge_types)
        
        volume_info[vol_tag] = (volume_info[vol_tag][0], periodic)

    for dim, vol_tag in volumes:
        print(f"\nVolume {vol_tag} classified as type {volume_info[vol_tag][0]}")
        surface_tags = gmsh.model.getBoundary([(dim, vol_tag)], 
                                        combined=False, oriented=False, recursive=False)
        for dim, surf_tag in surface_tags:    
            print(f"surface_info[{vol_tag}][{surf_tag}] ={surface_info[vol_tag][surf_tag]}")        
            edge_tags = gmsh.model.getBoundary([(dim, surf_tag)], 
                                        combined=False, oriented=False, recursive=False)
            for dim, edge_tag in edge_tags:    
                print(f"  edge_info[{vol_tag}][{surf_tag}][{edge_tag}] = {edge_info[vol_tag][surf_tag][edge_tag]}")
    
    return volume_info, surface_info, edge_info



def apply_transfinite_meshing10(volume_info, surface_info, edge_info, n_vertical, n_circumferential, n_radial, enable_periodic=False):

    surface_classes = {} 

    periodic_pairs = []
    
    edge2divisions = {}

    edges = gmsh.model.getEntities(1)
    surfaces = gmsh.model.getEntities(2)
    volumes = gmsh.model.getEntities(3)

    isPeriodic = True
    length = 0.0
    for dim, vol_tag in volumes:
        isPeriodic = isPeriodic and volume_info[vol_tag][1][0]
        length += volume_info[vol_tag][1][1] 
    length /= len(volumes)

    translation = [length, 0, 0 ]
    affine_transform = [
                1, 0, 0, translation[0],  # Row 1
                0, 1, 0, translation[1],  # Row 2
                0, 0, 1, translation[2],  # Row 3
                0, 0, 0, 1                # Row 4
            ]

    for dim, vol_tag in volumes:

        surface_tags = gmsh.model.getBoundary([(3, vol_tag)], 
                                        combined=False, oriented=False, recursive=False)

        inlet_tag = -1
        outlet_tag = -1
        for dim, surf_tag in surface_tags:
            surface_type = surface_info[vol_tag][surf_tag]
           
            if('inlet' in surface_type):
                surface_class_type = "inlet"
                inlet_tag = surf_tag
            elif('outlet' in surface_type):
                surface_class_type = "outlet"
                outlet_tag = surf_tag
            else:
                surface_class_type = "walls"


            if surface_class_type not in surface_classes:
                surface_classes[surface_class_type] = []
            surface_classes[surface_class_type].append(surf_tag)
            
            edge_classes = {} 
            edge_tags = gmsh.model.getBoundary([(2, surf_tag)], 
                                        combined=False, oriented=False, recursive=False)
            for dim, edge_tag in edge_tags:
                edge_type = edge_info[vol_tag][surf_tag][edge_tag]

                if edge_type not in edge_classes:
                    edge_classes[edge_type] = []
                edge_classes[edge_type].append(edge_tag)

                if edge_type == "along_edge":
                    #gmsh.model.mesh.setTransfiniteCurve(edge_tag, n_vertical)
                    edge2divisions[edge_tag] = n_vertical
                elif edge_type == "radial_edge":
                    #gmsh.model.mesh.setTransfiniteCurve(edge_tag, n_radial)
                    edge2divisions[edge_tag] = n_radial

                elif edge_type == "outer_circumferential_edge":
                    #gmsh.model.mesh.setTransfiniteCurve(edge_tag, n_circumferential)
                    edge2divisions[edge_tag] = n_circumferential

                elif edge_type == "inner_circumferential_edge":
                    #gmsh.model.mesh.setTransfiniteCurve(edge_tag, n_circumferential)
                    edge2divisions[edge_tag] = n_circumferential

            
            #for etype, edge_tags in edge_classes.items():
            ##    gmsh.model.addPhysicalGroup(1, edge_tags, name=etype)
        
            # gmsh.model.mesh.setTransfiniteSurface(surf_tag)
            # gmsh.model.mesh.setRecombine(2, surf_tag)

        if(enable_periodic and isPeriodic):
            gmsh.model.mesh.setPeriodic(2, [outlet_tag], [inlet_tag], affine_transform)
            print(f"    ✓ Applied periodic constraint: surface {outlet_tag} → {inlet_tag}")


    for edge_tag, divisions in edge2divisions.items():
        gmsh.model.mesh.setTransfiniteCurve(edge_tag, divisions)

            # for dim, edge_tag in edge_tags:
            #     edge_type = edge_info[vol_tag][surf_tag][edge_tag]

            #     if edge_type not in edge_classes:
            #         edge_classes[edge_type] = []
            #     edge_classes[edge_type].append(edge_tag)

            #     if edge_type == "along_edge":
            #         gmsh.model.mesh.setTransfiniteCurve(edge_tag, n_vertical)
            #     elif edge_type == "radial_edge":
            #         gmsh.model.mesh.setTransfiniteCurve(edge_tag, n_radial)
            #     elif edge_type == "outer_circumferential_edge":
            #         gmsh.model.mesh.setTransfiniteCurve(edge_tag, n_circumferential)
            #     elif edge_type == "inner_circumferential_edge":
            #         gmsh.model.mesh.setTransfiniteCurve(edge_tag, n_circumferential)

    for dim, surf_tag in surfaces:
        gmsh.model.mesh.setTransfiniteSurface(surf_tag)
        gmsh.model.mesh.setRecombine(2, surf_tag)
    
    for dim, vol_tag in volumes:
        gmsh.model.mesh.setTransfiniteVolume(vol_tag)

    for stype, surf_tags in surface_classes.items():
        print(f"Physical group surface: {len(surf_tags)} {stype} with tags {surf_tags}")
        gmsh.model.addPhysicalGroup(2, surf_tags, name=stype)

    gmsh.model.addPhysicalGroup(3, [vol_tag for dim, vol_tag in volumes], name="fluid")

    return isPeriodic, periodic_pairs


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

def define_physical_groups1(volume_info, surface_info, edge_info, enable_periodic=False):
    print("\nDefining physical groups for boundaries.")    


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

    # Initialize GMSH
    gmsh.initialize()
    gmsh.model.add("cyl_hex")

    gmsh.option.setNumber("General.Terminal", 1)
    
    # Import the STEP file
    print(f"Importing {step_file}...")
    #gmsh.merge(step_file)
    
# Import the STEP and build OCC shapes
    gmsh.model.occ.importShapes(step_file)
    gmsh.model.occ.synchronize()

    # 1) Remove duplicated geometric entities (same location, different tags)
    #    This helps avoid overlapping faces/edges coming from independent solids.
    try:
        gmsh.model.occ.removeAllDuplicates()
    except Exception:
        # Older API name
        gmsh.model.occ.removeDuplicateEntities()
    # Synchronize to get the model
    gmsh.model.occ.synchronize()

    vols = gmsh.model.getEntities(3)
    if vols:
        frags, _ = gmsh.model.occ.fragment(vols, [])
        gmsh.model.occ.synchronize()
    else:
        print("[WARN] No volumes found after STEP import.")

    # Increase the tolerance for matching points significantly
    # The geometry has some numerical precision issues from STEP import
    gmsh.option.setNumber("Geometry.Tolerance", 1e-5)
    gmsh.option.setNumber("Geometry.ToleranceBoolean", 1e-5)
    gmsh.option.setNumber("Mesh.ToleranceInitialDelaunay", 1e-5)
    

    volume_info, surface_info, edge_info =label_entities(epsilon=1.0e-8)

    inlet_plane_surfaces = []
    for vol_tag in surface_info.keys():
        for surf_tag in surface_info[vol_tag].keys():
            if surface_info[vol_tag][surf_tag] == "inlet_plane":
                inlet_plane_surfaces.append(surf_tag)
    gmsh.model.addPhysicalGroup(2, inlet_plane_surfaces, name="inlet")
    print("inlet_plane_surfaces:", inlet_plane_surfaces)

    outlet_plane_surfaces = []
    for vol_tag in surface_info.keys():
        for surf_tag in surface_info[vol_tag].keys():
            if surface_info[vol_tag][surf_tag] == "outlet_plane":
                outlet_plane_surfaces.append(surf_tag)
    gmsh.model.addPhysicalGroup(2, outlet_plane_surfaces, name="outlet")

    wall_plane_surfaces = []
    for vol_tag in surface_info.keys():
        for surf_tag in surface_info[vol_tag].keys():
            if surface_info[vol_tag][surf_tag] == "outer_circumferential_wall" :
                wall_plane_surfaces.append(surf_tag)
    gmsh.model.addPhysicalGroup(2, wall_plane_surfaces, name="wall")

    gmsh.model.addPhysicalGroup(3, [vol_tag for vol_tag in volume_info.keys()], name="fluid")


    # Apply transfinite meshing
    #apply_transfinite_meshing(n_vertical, n_circumferential, n_radial, surface_info)
    isPeriodic, periodic_pairs = apply_transfinite_meshing1(volume_info, surface_info, edge_info,n_vertical, n_circumferential, n_radial, enable_periodic=periodic)
    
    # Define physical groups and check periodicity
    #print("\nDefining physical groups:")
    #periodic_pairs = define_physical_groups1(volume_info, surface_info, edge_info, enable_periodic=periodic)
    
    # Apply periodic constraints BEFORE mesh generation
    if periodic and isPeriodic and periodic_pairs:
        apply_periodic_constraints(periodic_pairs)
    
    # Generate mesh
    print("\nGenerating mesh...")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    gmsh.model.mesh.generate(3)
    
    # Set mesh order to 2 for curved elements
    print("Creating high-order curved mesh...")
    #gmsh.model.mesh.setOrder(2)
    
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
