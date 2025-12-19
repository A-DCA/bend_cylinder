#!/usr/bin/env python3
"""
CadQuery script to view STEP files in cq-editor
Usage: Open this file in cq-editor
"""

import cadquery as cq

# Import the STEP file
result = cq.importers.importStep("swept_inner_outer_cylinder_curved.step")

# For cq-editor to display it, assign to 'show_object' or just use variable name
show_object(result)
