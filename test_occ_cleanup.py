#!/usr/bin/env python3
"""
Hex mesh for a multi-solid cylinder imported from STEP, using transfinite meshing.
- Unify duplicate OCC entities and fragment volumes to share interfaces.
- Apply transfinite constraints on curves/surfaces/volumes (adjust counts below).
- Force first-order elements and MSH v2.2 ASCII for gmshToFoam.

Usage:
    python gmsh_fix_hex_cylinder.py path/to/model.step

Notes:
- You need gmsh>=4.9 installed with the Python API.
- The STEP should contain the solids you described (inner core + 4 outer blocks x axial slices).
- Tweak NR, NT, NZ to your desired divisions.
"""
import sys
import gmsh

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


def apply_transfinite_meshing(volume_info, surface_info, edge_info, n_vertical, n_circumferential, n_radial, enable_periodic=False):

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

def addPhysicalGroups(volume_info, surface_info, edge_info):
    surface_classes = {"inlet": [], "outlet": [], "wall": []}
    for vol_tag in surface_info.keys():
        for surf_tag in surface_info[vol_tag].keys():
            if surface_info[vol_tag][surf_tag] == "inlet_plane":
                surface_classes["inlet"].append(surf_tag)
            if surface_info[vol_tag][surf_tag] == "outlet_plane":
                surface_classes["outlet"].append(surf_tag)
            if surface_info[vol_tag][surf_tag] == "outer_circumferential_wall" :
                surface_classes["wall"].append(surf_tag)

    gmsh.model.addPhysicalGroup(2, surface_classes["wall"], name="walls")    
    gmsh.model.addPhysicalGroup(2, surface_classes["inlet"], name="inlet")
    gmsh.model.addPhysicalGroup(2, surface_classes["outlet"], name="outlet")

    gmsh.model.addPhysicalGroup(3, [vol_tag for vol_tag in volume_info.keys()], name="internal")

def setPeriodicBs(volume_info, surface_info, edge_info, enable_periodic):
    isPeriodic = True
    length = 0.0
    for vol_tag in volume_info.keys():
        isPeriodic = isPeriodic and volume_info[vol_tag][1][0]
        length += volume_info[vol_tag][1][1] 
    length /= len(volume_info.keys())

    translation = [length, 0, 0 ]
    affine_transform = [
                1, 0, 0, translation[0],  # Row 1
                0, 1, 0, translation[1],  # Row 2
                0, 0, 1, translation[2],  # Row 3
                0, 0, 0, 1                # Row 4
            ]
    if(enable_periodic and isPeriodic):
        for vol_tag in volume_info.keys():
            inlet_tag = -1
            outlet_tag = -1
            surface_tags = gmsh.model.getBoundary([(3, vol_tag)], 
                                        combined=False, oriented=False, recursive=False)
            for dim, surf_tag in surface_tags:
                surface_type = surface_info[vol_tag][surf_tag]
               
                if('inlet' in surface_type):
                    inlet_tag = surf_tag
                elif('outlet' in surface_type):
                    outlet_tag = surf_tag

            gmsh.model.mesh.setPeriodic(2, [outlet_tag], [inlet_tag], affine_transform)
            print(f"    ✓ Applied periodic constraint: surface {outlet_tag} → {inlet_tag}")

def setTransfiniteMeshing(volume_info, surface_info, edge_info, NL, NT, NR):
    for vol_tag in surface_info.keys():
        for surf_tag in surface_info[vol_tag].keys():
            for edge_tag in edge_info[vol_tag][surf_tag].keys():
                edge_type = edge_info[vol_tag][surf_tag][edge_tag]
                ndividions = 1
                if edge_type == "along_edge":
                    ndividions = NL
                elif edge_type == "radial_edge":
                    ndividions = NR
                elif edge_type == "outer_circumferential_edge":
                    ndividions = NT
                elif edge_type == "inner_circumferential_edge":
                    ndividions = NT
                #edge
                gmsh.model.mesh.setTransfiniteCurve(edge_tag, ndividions)
            #surface    
            gmsh.model.mesh.setTransfiniteSurface(surf_tag)
            # Recombine to quads on surfaces (required precursor for hexes)
            gmsh.model.mesh.setRecombine(2, surf_tag)
        #volume
        gmsh.model.mesh.setTransfiniteVolume(vol_tag)

    # # 3) Apply transfinite constraints
    # #    We set a uniform number of points on curves and mark surfaces/volumes.
    # #    If needed, refine only selected entities by editing below.
    # for dim, tag in gmsh.model.getEntities(1):  # curves
    #     # Heuristic: axial vs circumferential vs radial classification is tricky
    #     # on imported STEP; use a single count to keep topology consistent.
    #     #gmsh.model.mesh.setTransfiniteCurve(tag, max(NR, NT, NZ))
    #     gmsh.model.mesh.setTransfiniteCurve(tag, edge2divisions[tag])


    # for dim, tag in gmsh.model.getEntities(2):  # surfaces
    #     gmsh.model.mesh.setTransfiniteSurface(tag)
    #     # Recombine to quads on surfaces (required precursor for hexes)
    #     gmsh.model.mesh.setRecombine(dim, tag)

    # for dim, tag in gmsh.model.getEntities(3):  # volumes
    #     gmsh.model.mesh.setTransfiniteVolume(tag)

# ---- User-tunable transfinite counts ----
NR = 10   # radial divisions per block
NT = 10  # circumferential divisions per quadrant (outer blocks)
NL = 10  # axial divisions across the cylinder length - longanituional 
enable_periodic = True
# -----------------------------------------

step_path = sys.argv[1] if len(sys.argv) > 1 else "swept_inner_outer_cylinder_straight.step"

gmsh.initialize()
try:
    gmsh.model.add("cyl_hex")

    # Import the STEP and build OCC shapes
    gmsh.model.occ.importShapes(step_path)
    gmsh.model.occ.synchronize()

    # 1) Remove duplicated geometric entities (same location, different tags)
    #    This helps avoid overlapping faces/edges coming from independent solids.
    try:
        gmsh.model.occ.removeAllDuplicates()
    except Exception:
        # Older API name
        gmsh.model.occ.removeDuplicateEntities()
    gmsh.model.occ.synchronize()

    # 2) Fragment all volumes together so that interfaces are shared
    vols = gmsh.model.getEntities(3)
    if vols:
        frags, _ = gmsh.model.occ.fragment(vols, [])
        gmsh.model.occ.synchronize()
    else:
        print("[WARN] No volumes found after STEP import.")

    print(f"{len(gmsh.model.getEntities(3))} volumes to be processed.")
    print(f"{len(gmsh.model.getEntities(2))} surfaces to be processed.")
    print(f"{len(gmsh.model.getEntities(1))} edges to be processed.")
    print(f"{len(gmsh.model.getEntities(0))} vertices to be processed.")

    volume_info, surface_info, edge_info = label_entities()
    
    addPhysicalGroups(volume_info, surface_info, edge_info)

    # surface_classes = {"inlet": [], "outlet": [], "wall": []}
    # for vol_tag in surface_info.keys():
    #     for surf_tag in surface_info[vol_tag].keys():
    #         if surface_info[vol_tag][surf_tag] == "inlet_plane":
    #             surface_classes["inlet"].append(surf_tag)
    #         if surface_info[vol_tag][surf_tag] == "outlet_plane":
    #             surface_classes["outlet"].append(surf_tag)
    #         if surface_info[vol_tag][surf_tag] == "outer_circumferential_wall" :
    #             surface_classes["wall"].append(surf_tag)

            
    # gmsh.model.addPhysicalGroup(2, surface_classes["inlet"], name="inlet")
    # gmsh.model.addPhysicalGroup(2, surface_classes["outlet"], name="outlet")
    # gmsh.model.addPhysicalGroup(2, surface_classes["wall"], name="wall")
    # gmsh.model.addPhysicalGroup(3, [vol_tag for vol_tag in volume_info.keys()], name="fluid")

    setPeriodicBs(volume_info, surface_info, edge_info, enable_periodic)
    # isPeriodic = True
    # length = 0.0
    # for vol_tag in volume_info.keys():
    #     isPeriodic = isPeriodic and volume_info[vol_tag][1][0]
    #     length += volume_info[vol_tag][1][1] 
    # length /= len(volume_info.keys())

    # translation = [length, 0, 0 ]
    # affine_transform = [
    #             1, 0, 0, translation[0],  # Row 1
    #             0, 1, 0, translation[1],  # Row 2
    #             0, 0, 1, translation[2],  # Row 3
    #             0, 0, 0, 1                # Row 4
    #         ]
    # if(enable_periodic and isPeriodic):
    #     for vol_tag in volume_info.keys():
    #         inlet_tag = -1
    #         outlet_tag = -1
    #         surface_tags = gmsh.model.getBoundary([(3, vol_tag)], 
    #                                     combined=False, oriented=False, recursive=False)
    #         for dim, surf_tag in surface_tags:
    #             surface_type = surface_info[vol_tag][surf_tag]
               
    #             if('inlet' in surface_type):
    #                 inlet_tag = surf_tag
    #             elif('outlet' in surface_type):
    #                 outlet_tag = surf_tag

    #         gmsh.model.mesh.setPeriodic(2, [outlet_tag], [inlet_tag], affine_transform)
    #         print(f"    ✓ Applied periodic constraint: surface {outlet_tag} → {inlet_tag}")

    setTransfiniteMeshing(volume_info, surface_info, edge_info, NL, NT, NR)    
    # edge2divisions  = {}
    # for vol_tag in surface_info.keys():
    #     for surf_tag in surface_info[vol_tag].keys():
    #         for edge_tag in edge_info[vol_tag][surf_tag].keys():
    #             edge_type = edge_info[vol_tag][surf_tag][edge_tag]
    #             if edge_type == "along_edge":
    #                 edge2divisions[edge_tag] = NL
    #             elif edge_type == "radial_edge":
    #                 edge2divisions[edge_tag] = NR
    #             elif edge_type == "outer_circumferential_edge":
    #                 edge2divisions[edge_tag] = NT
    #             elif edge_type == "inner_circumferential_edge":
    #                 edge2divisions[edge_tag] = NT
    # # 3) Apply transfinite constraints
    # #    We set a uniform number of points on curves and mark surfaces/volumes.
    # #    If needed, refine only selected entities by editing below.
    # for dim, tag in gmsh.model.getEntities(1):  # curves
    #     # Heuristic: axial vs circumferential vs radial classification is tricky
    #     # on imported STEP; use a single count to keep topology consistent.
    #     #gmsh.model.mesh.setTransfiniteCurve(tag, max(NR, NT, NZ))
    #     gmsh.model.mesh.setTransfiniteCurve(tag, edge2divisions[tag])


    # for dim, tag in gmsh.model.getEntities(2):  # surfaces
    #     gmsh.model.mesh.setTransfiniteSurface(tag)
    #     # Recombine to quads on surfaces (required precursor for hexes)
    #     gmsh.model.mesh.setRecombine(dim, tag)

    # for dim, tag in gmsh.model.getEntities(3):  # volumes
    #     gmsh.model.mesh.setTransfiniteVolume(tag)



    # 4) Mesh options for OpenFOAM
    gmsh.option.setNumber("Mesh.ElementOrder", 1)          # first-order only
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)      # MSH v2.2
    gmsh.option.setNumber("Mesh.Binary", 0)                # ASCII
    gmsh.option.setNumber("Mesh.RecombineAll", 1)          # global recombination
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1) # Simple algorithm

    # 5) Generate 3D mesh and write
    gmsh.model.mesh.generate(3)

    outname = "cylinder_hex_v2.msh"
    gmsh.write(outname)
    print(f"Wrote: {outname}")

finally:
    gmsh.finalize()
