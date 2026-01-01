#!/usr/bin/env python3
"""
OpenFOAM Periodic Tube Case Setup Script
Automates the workflow from STEP file to running simulation with periodic boundaries
"""

import os
import re
import sys
import subprocess
import shutil
from pathlib import Path
import json

from matplotlib import patches

def get_boundary_patches(boundary_file):
 # Read boundary file (adjust the path if needed)
    text = Path(boundary_file).read_text(encoding='utf-8')

    # Remove // comments
    clean = re.sub(r"//.*", "", text)

    # Capture patch blocks: name { ... }
    pattern = re.compile(r"\b([A-Za-z0-9_]+)\s*\{([^}]*)\}", re.MULTILINE | re.DOTALL)
    patches = []
    for name, body in pattern.findall(clean):
        entries = {}
        for line in body.splitlines():
            line = line.strip()
            if not line:
                continue
            # Lines of the form: key value; or key (vector); or key string;
            m = re.match(r"([A-Za-z0-9_]+)\s+(.+?);?$", line)
            if not m:
                continue
            key, val = m.groups()
            val = val.strip()

            # Vector like (0.1 0 0)
            if re.match(r"^\([^)]+\)$", val):
                parts = val.strip('()').split()
                try:
                    entries[key] = [float(x) for x in parts]
                except ValueError:
                    entries[key] = parts
            # Integer
            elif re.match(r"^[0-9]+$", val):
                entries[key] = int(val)
            # Float (supports scientific notation)
            elif re.match(r"^[0-9]+(\.[0-9]+)?([Ee][\+\-]?[0-9]+)?$", val):
                entries[key] = float(val)
            else:
                # Raw string (e.g. 'translational', 'outlet', 'ascii', etc.)
                entries[key] = val

        patches.append({"name": name, **entries})

    result = {"nPatches": len(patches), "patches": patches}
    return result

class PeriodicTubeCaseSetup:
    """Setup and run OpenFOAM periodic tube flow simulation"""
    
    def __init__(self, case_dir, step_file, gmsh_file, mesh_params, flow_params):
        """
        Initialize case setup
        
        Args:
            case_dir: Path to case directory
            step_file: Path to STEP geometry file
            params: Dictionary with simulation parameters
        """

        self.step_file =""
        self.gmsh_file =""
        
        self.case_dir = Path(case_dir)
        # Create case directory if it doesn't exist
        self.case_dir.mkdir(parents=True, exist_ok=True)

        if( step_file is None and gmsh_file is None):
            raise ValueError("Either step_file or gmsh_file must be provided.")
        elif( step_file is None):
            self.gmsh_file = Path(gmsh_file)
            shutil.copy(self.gmsh_file, self.case_dir)
        elif(gmsh_file is None):
            self.step_file = Path(step_file)
        else:
            self.step_file = Path(step_file)
            self.gmsh_file = Path(gmsh_file)
            shutil.copy(self.gmsh_file, self.case_dir)

        self.mesh_params = mesh_params
        self.flow_params = flow_params
        
        
    def run_command(self, cmd, description, cwd=None):
        """Execute shell command with error handling"""
        print(f"\n{'='*70}")
        print(f"STEP: {description}")
        print(f"{'='*70}")
        print(f"Command: {cmd}")
        
        working_dir = cwd or self.case_dir
        result = subprocess.run(
            cmd, 
            shell=True, 
            cwd=working_dir,
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            print(f"ERROR: {description} failed!")
            print(f"stderr: {result.stderr}")
            return False
        
        print(f"SUCCESS: {description} completed")
        if result.stdout:
            print(f"Output: {result.stdout[-500:]}")  # Last 500 chars
        return True
    
    
    
    def convert_mesh(self, mesh_file):
        """Convert gmsh mesh to OpenFOAM format"""
        
        cmd = f"gmshToFoam -keepOrientation {mesh_file.name}"
        if not self.run_command(cmd, "Convert mesh to OpenFOAM format"):
            return False
        
        return True
    
    def identify_boundaries(self):
        """Identify and name boundary patches"""
        
        cmd = "autoPatch 30 -overwrite"
        if not self.run_command(cmd, "Identify boundary patches"):
            return False
        
        return True
    def setup_cyclicAMI_boundaries(self):
        boundary_file = self.case_dir / "constant/polyMesh/boundary"
        patches_info = get_boundary_patches(boundary_file)    
        print(f"Found {patches_info['nPatches']} patches in boundary file.")

        L = self.mesh_params['length']
        translate = self.mesh_params['translate']

        all_patches = []
        patch_dict = {p['name']: p for p in patches_info['patches']}

        required_patches = {
            'inlet': 'cyclicAMI',
            'outlet': 'cyclicAMI',
            'walls': 'wall'
        }
        for name, expected_type in required_patches.items():
            patch = patch_dict.get(name)
            if patch is None:
                raise ValueError(f"'{name}' patch not found.")
            for key in ['nFaces', 'startFace']:
                if key not in patch:
                    raise ValueError(f"'{key}' not found in {name} patch.")
            content = ""
            if( name == 'inlet'):
                l = 1
                frm = 'inlet'
                to = 'outlet'
            elif( name == 'outlet'):
                l = -1
                to = 'inlet'
                frm = 'outlet'
            if( name != 'walls'):
                content = f'''{frm}
                    {{
                        type            {expected_type};
                        matchTolerance  0.01;
                        transform       translational;
                        neighbourPatch  {to};
                        separationVector ({l*translate[0]} {l*translate[1]} {l*translate[2]});
                        nFaces          {patch['nFaces']};
                        startFace       {patch['startFace']};
                    }}'''
            else:
                content = f'''{name}
                    {{
                        type            {expected_type};
                        inGroups        1(wall);
                        nFaces          {patch['nFaces']};
                        startFace       {patch['startFace']};
                    }}'''
            
            all_patches.append((name, patch, content))

# Sort by startFace
        all_patches.sort(key=lambda x: x[1]['startFace'])
        
        # Build boundary file content
        patches_str = '\n    '.join([p[2] for p in all_patches])
        
        # Write new boundary file
        header = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2506                                  |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       polyBoundaryMesh;
    location    "constant/polyMesh";
    object      boundary;
    author      "jma";
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

3
(
    {patches_str}
)

// ************************************************************************* //
"""
        
        with open(boundary_file, 'w') as f:
            f.write(header)
        
        print("Modified boundary file with cyclicAMI patches")
        return True
                

    def setup_cyclicAMI_boundaries0(self):
        """Modify boundary file to set up cyclicAMI patches and rename patches"""
        
        boundary_file = self.case_dir / "constant/polyMesh/boundary"

        patches_info = get_boundary_patches(boundary_file)
        print(f"Found {patches_info['nPatches']} patches in boundary file.")
        # Backup original
        shutil.copy(boundary_file, str(boundary_file) + '.backup')
        
        # Use changeDictionary or manual rewrite with better parsing
        # Read the backup to parse properly
        with open(str(boundary_file) + '.backup', 'r') as f:
            lines = f.readlines()
        
        L = self.params['length']
        
        # Parse patches manually - look for pattern: patchName { ... nFaces ... startFace ... }
        patches_data = []
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            # Check if next line has opening brace (indicates patch definition)
            if (i + 1 < len(lines) and '{' in lines[i+1] and 
                line and not line.startswith('//') and not line.startswith('/*') and
                'FoamFile' not in line and '(' not in line and ')' not in line):
                
                patch_name = line
                patch_info = {'name': patch_name}
                i += 2  # Skip name and opening brace
                
                # Parse patch contents until closing brace
                while i < len(lines) and '}' not in lines[i]:
                    if 'nFaces' in lines[i]:
                        patch_info['nFaces'] = int(lines[i].split()[1].rstrip(';'))
                    if 'startFace' in lines[i]:
                        patch_info['startFace'] = int(lines[i].split()[1].rstrip(';'))
                    i += 1
                
                if 'nFaces' in patch_info and 'startFace' in patch_info and patch_info['nFaces'] > 0:
                    patches_data.append(patch_info)
                    print(f"Found patch: {patch_info['name']} with {patch_info['nFaces']} faces")
            i += 1
        
        if len(patches_data) < 3:
            print(f"ERROR: Found only {len(patches_data)} patches")
            return False
        
        """This sorting below results mismatched patches sometimes. using order found instead."""
        # Sort by nFaces to identify patches
        #sorted_patches = sorted(patches_data, key=lambda x: x['nFaces'])
        
        sorted_patches = patches_data
        inlet_info = sorted_patches[0]
        outlet_info = sorted_patches[1]
        walls_info = sorted_patches[2]
        
        print(f"Identified: inlet={inlet_info['name']} ({inlet_info['nFaces']} faces), "
              f"outlet={outlet_info['name']} ({outlet_info['nFaces']} faces), "
              f"walls={walls_info['name']} ({walls_info['nFaces']} faces)")
        
        # Sort patches by startFace to maintain correct order in boundary file
        all_patches = [
            ('inlet', inlet_info, f'''inlet
    {{
        type            cyclicAMI;
        matchTolerance  0.01;
        transform       translational;
        neighbourPatch  outlet;
        separationVector ({L} 0 0);
        nFaces          {inlet_info['nFaces']};
        startFace       {inlet_info['startFace']};
    }}'''),
            ('walls', walls_info, f'''walls
    {{
        type            wall;
        inGroups        1(wall);
        nFaces          {walls_info['nFaces']};
        startFace       {walls_info['startFace']};
    }}'''),
            ('outlet', outlet_info, f'''outlet
    {{
        type            cyclicAMI;
        matchTolerance  0.01;
        transform       translational;
        neighbourPatch  inlet;
        separationVector ({-L} 0 0);
        nFaces          {outlet_info['nFaces']};
        startFace       {outlet_info['startFace']};
    }}''')
        ]
        
        # Sort by startFace
        all_patches.sort(key=lambda x: x[1]['startFace'])
        
        # Build boundary file content
        patches_str = '\n    '.join([p[2] for p in all_patches])
        
        # Write new boundary file
        header = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2506                                  |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       polyBoundaryMesh;
    location    "constant/polyMesh";
    object      boundary;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

3
(
    {patches_str}
)

// ************************************************************************* //
"""
        
        with open(boundary_file, 'w') as f:
            f.write(header)
        
        print("Modified boundary file with cyclicAMI patches")
        return True
    
    def create_initial_conditions(self):
        """Create 0/ directory with initial conditions"""
        
        # Delete existing 0 directory if present
        zero_dir = self.case_dir / "0"
        if zero_dir.exists():
            shutil.rmtree(zero_dir)
        zero_dir.mkdir(exist_ok=True)
        
        # Velocity field (0/U)
        U_file = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{{
    walls
    {{
        type            noSlip;
    }}

    inlet
    {{
        type            cyclicAMI;
        value           uniform (0 0 0);
    }}

    outlet
    {{
        type            cyclicAMI;
        value           uniform (0 0 0);
    }}
}}

// ************************************************************************* //
"""
        
        # Pressure field (0/p)
        p_file = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{{
    walls
    {{
        type            zeroGradient;
    }}

    inlet
    {{
        type            cyclicAMI;
        value           uniform 0;
    }}

    outlet
    {{
        type            cyclicAMI;
        value           uniform 0;
    }}
}}

// ************************************************************************* //
"""
        
        with open(zero_dir / "U", 'w') as f:
            f.write(U_file)
        
        with open(zero_dir / "p", 'w') as f:
            f.write(p_file)
        
        print("Created initial condition files (0/U, 0/p)")
        return True
    
    def create_transport_properties(self):
        """Create constant/transportProperties"""
        
        const_dir = self.case_dir / "constant"
        const_dir.mkdir(exist_ok=True)
        
        nu = self.flow_params['kinematic_viscosity']
        
        transport_file = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

transportModel  Newtonian;

nu              {nu};

// ************************************************************************* //
"""
        
        with open(const_dir / "transportProperties", 'w') as f:
            f.write(transport_file)
        
        print("Created constant/transportProperties")
        return True
    
    def create_turbulence_properties(self):
        """Create constant/turbulenceProperties"""
        
        const_dir = self.case_dir / "constant"
        
        turb_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType  laminar;

// ************************************************************************* //
"""
        
        with open(const_dir / "turbulenceProperties", 'w') as f:
            f.write(turb_file)
        
        print("Created constant/turbulenceProperties")
        return True
    
    def create_fvOptions(self):
        """Create constant/fvOptions for pressure gradient source"""
        
        const_dir = self.case_dir / "constant"
        
        dp_dx = self.flow_params['pressure_gradient']
        rho = self.flow_params.get('density', 1000)  # Default water density
        
        # Convert pressure gradient to acceleration: a = (dp/dx) / rho
        accel = dp_dx / rho
        
        fvOptions_file = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvOptions;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

pressureGradient
{{
    type            vectorSemiImplicitSource;
    
    selectionMode   all;
    volumeMode      specific;
    
    sources
    {{
        U           (({accel} 0 0) 0);
    }}
}}

// ************************************************************************* //
"""
        
        with open(const_dir / "fvOptions", 'w') as f:
            f.write(fvOptions_file)
        
        print(f"Created constant/fvOptions (pressure gradient: {dp_dx} Pa/m)")
        return True
    
    def create_control_dict(self):
        """Create system/controlDict"""
        
        sys_dir = self.case_dir / "system"
        sys_dir.mkdir(exist_ok=True)
        
        end_time = self.flow_params.get('end_time', 100)
        write_interval = self.flow_params.get('write_interval', 100)
        
        control_file = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     simpleFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         {end_time};

deltaT          1;

writeControl    timeStep;

writeInterval   {write_interval};

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

// ************************************************************************* //
"""
        
        with open(sys_dir / "controlDict", 'w') as f:
            f.write(control_file)
        
        print("Created system/controlDict")
        return True
    
    def create_fvSchemes(self):
        """Create system/fvSchemes"""
        
        sys_dir = self.case_dir / "system"
        
        schemes_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

// ************************************************************************* //
"""
        
        with open(sys_dir / "fvSchemes", 'w') as f:
            f.write(schemes_file)
        
        print("Created system/fvSchemes")
        return True
    
    def create_fvSolution(self):
        """Create system/fvSolution"""
        
        sys_dir = self.case_dir / "system"
        
        solution_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-08;
        relTol          0.01;
        smoother        GaussSeidel;
    }

    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-08;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 2;
    consistent      yes;
    
    pRefCell        0;
    pRefValue       0;
    
    residualControl
    {
        p               1e-5;
        U               1e-5;
    }
}

relaxationFactors
{
    fields
    {
        p               0.2;
    }
    equations
    {
        U               0.5;
    }
}

// ************************************************************************* //
"""
        
        with open(sys_dir / "fvSolution", 'w') as f:
            f.write(solution_file)
        
        print("Created system/fvSolution")
        return True
    
    def check_mesh(self):
        """Run checkMesh to verify mesh quality"""
        
        cmd = "checkMesh"
        return self.run_command(cmd, "Check mesh quality")
    
    def run_simulation(self):
        """Run simpleFoam solver"""
        
        cmd = "simpleFoam 2>&1 | tee log.simpleFoam"
        return self.run_command(cmd, "Run simpleFoam simulation")
    
    def setup_and_run(self):
        """Execute complete workflow"""
        
        print("\n" + "="*70)
        print("OpenFOAM Periodic Tube Case Setup")
        print("="*70)
        print(f"Case directory: {self.case_dir}")
        print(f"Gmsh file: {self.gmsh_file}")
        print(f"Parameters: {self.mesh_params}")
        print(f"Parameters: {self.flow_params}")
        
        # Step 1: Create system files first (needed by OpenFOAM commands)
        print("\n[1/11] Creating system files...")
        if not self.create_control_dict():
            return False
        if not self.create_fvSchemes():
            return False
        if not self.create_fvSolution():
            return False
        
        # Step 2: Generate mesh
        print("\n[2/11] Generating mesh with gmsh...")
        if self.gmsh_file.exists():
            print(f"Gmsh file {self.gmsh_file} exists!")
            mesh_file = self.gmsh_file
        else:
            print(f"no Gmsh mesh file {self.gmsh_file} found.")
            return False
        
        # Step 3: Convert mesh
        print("\n[3/11] Converting mesh to OpenFOAM format...")
        if not self.convert_mesh(mesh_file):
            return False
        
        # Step 4: Identify boundaries
        #print("\n[4/11] Identifying boundary patches...")
        #if not self.identify_boundaries():
        #    return False
        
        # Step 5: Setup cyclicAMI boundaries
        print("\n[5/11] Setting up cyclicAMI boundaries...")
        if not self.setup_cyclicAMI_boundaries():
            return False
        
        # Step 6: Create initial conditions
        print("\n[6/11] Creating initial conditions...")
        if not self.create_initial_conditions():
            return False
        
        # Step 7: Create transport properties
        print("\n[7/11] Creating transport properties...")
        if not self.create_transport_properties():
            return False
        
        # Step 8: Create turbulence properties
        print("\n[8/11] Creating turbulence properties...")
        if not self.create_turbulence_properties():
            return False
        
        # Step 9: Create fvOptions
        print("\n[9/11] Creating fvOptions (pressure gradient source)...")
        if not self.create_fvOptions():
            return False
        
        # Step 10: Check mesh
        print("\n[10/11] Checking mesh quality...")
        if not self.check_mesh():
            print("WARNING: Mesh check failed, but continuing...")
        
        # Step 11: Run simulation
        print("\n[11/11] Running simulation...")
        if not self.run_simulation():
            return False
        
        print("\n" + "="*70)
        print("SETUP AND SIMULATION COMPLETE!")
        print("="*70)
        print(f"\nTo visualize results, run: paraFoam -case {self.case_dir}")
        print("\nTo compute physical pressure in ParaView, use Calculator filter:")
        print(f"  p_physical = p - {self.flow_params['pressure_gradient']} * {self.mesh_params['length']} ")
        
        return True


def main():
    """Main execution function"""
    if len(sys.argv) < 3:
        print("Usage: python setup_periodic_tube.py mesh_config_file mesh_file")
        sys.exit(1)

    mesh_config_file = sys.argv[1]
    mesh_file = sys.argv[2]
    if not mesh_file.lower().endswith('.msh') and not Path.exists(mesh_file):
        print("Error: The first argument must be an existing Gmsh mesh file.")
        sys.exit(1)

    #check mesh periodic in gmsh
    with open(mesh_file, 'r', encoding='utf-8', errors='ignore') as f:
        mesh_contents = f.read()
    if '$Periodic' not in mesh_contents:
        print("Warning: The mesh file does not contain any $Periodic section. Periodic boundaries may not be defined in the mesh.")
        sys.exit(1)
    else:
        print("Mesh file contains $Periodic section.")
    

    # Example parameters
    flow_params = {
        'kinematic_viscosity': 1.0e-3,    # Kinematic viscosity [m²/s] (water)
        'density': 1000,                  # Density [kg/m³] (water)
        'pressure_gradient': 10000,        # Pressure gradient [Pa/m]
        'end_time': 200,                 # Number of iterations
        'write_interval': 200,           # Write interval
    }
    

    with open(mesh_config_file, 'r') as f:
        mesh_params = json.load(f)

    case_name = os.path.basename(mesh_file).split('.')[0]
    case_dir = Path.cwd()/ case_name
    
    # Create and run setup
    setup = PeriodicTubeCaseSetup(case_dir, None, mesh_file, mesh_params, flow_params)
    success = setup.setup_and_run()
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
