const char * objcurves_py = R"MULTILINE(
# Blender Python API script used to convert OBJ lines into curves
# to give them material with color.
# In Blender 4.0, OBJ lines already appear to be loaded as curves.
# Either that, or all the objects imported were still selected
# when conversion was carried out.
# The JSON file name placeholder in this file is replaced with
# the actual path to the JSON file produced by the Netmorph OBJ
# generator, and the resulting temporary Python script is given
# to Blender.

import bpy
from mathutils import *
import json

def get_object_info()->tuple:
	# Obtain necessary information about the OBJ file to properly
	# run the conversion through Blender.
	with open('OBJDATA_JSON', 'r') as f:
		obj_data = json.load(f)

	filepath = obj_data['obj_path']
	axons_list = obj_data['axons']
	dendrites_list = obj_data['dendrites']
	somas_list = obj_data['somas']
	synapses_list = obj_data['synapses']
	blendpath = obj_data['blend_path']
	axonbevdepth = obj_data['axonbevdepth']
	dendritebevdepth = obj_data['dendritebevdepth']
	return filepath, axons_list, dendrites_list, somas_list, synapses_list, blendpath, axonbevdepth, dendritebevdepth

def get_blender_objects(names_list:list)->list:
	blender_objects = []
	for name in names_list:
		blender_objects.append( bpy.data.objects[name] )
	return blender_objects

def set_bevel_depths(objects_list:list, depth:float):
	for curve_object in objects_list:
		curve_object.data.bevel_depth = depth

def delete_default_cube():
	# Delete the default Cube object
	# 1. Deselect all
	bpy.ops.object.select_all(action='DESELECT')
	# 3. Select the Cube
	bpy.data.objects['Cube'].select_set(True)
	# 4. Delete it    
	bpy.ops.object.delete()

def add_material_for_object(blender_object, object_name:str, color:tuple):
	#create new material with name of object
    new_mat = bpy.data.materials.new(object_name)
    #add new material to object
    blender_object.data.materials.clear()
    blender_object.data.materials.append(new_mat)
    new_mat.diffuse_color = color

filepath, axons_list, dendrites_list, somas_list, synapses_list, blendpath, axonbevdepth, dendritebevdepth = get_object_info()

# Note: After importing, all objects are automatically selected
# so that converting to curves converts all of the lines to
# curves without having to select them individually to do so.
bpy.ops.wm.obj_import(filepath=filepath)

# Convert all selected to curves
bpy.ops.object.convert(target='CURVE')

axon_objects = get_blender_objects(axons_list)
dendrite_objects = get_blender_objects(dendrites_list)
soma_objects = get_blender_objects(somas_list)
synapses_objects = get_blender_objects(synapses_list)

set_bevel_depths(axon_objects, axonbevdepth)
set_bevel_depths(dendrite_objects, dendritebevdepth)

for i in range(len(axons_list)):
	axon_name = axons_list[i]
	axon_object = axon_objects[i]
	add_material_for_object(axon_object, axon_name, (0.0, 1.0, 0.0, 1.0))

for i in range(len(dendrites_list)):
	dendrite_name = dendrites_list[i]
	dendrite_object = dendrite_objects[i]
	add_material_for_object(dendrite_object, axon_name, (0.0, 0.0, 1.0, 1.0))

for i in range(len(synapses_list)):
	synapse_name = synapses_list[i]
	synapse_object = synapses_objects[i]
	add_material_for_object(synapse_object, synapse_name, (1.0, 0.0, 0.0, 1.0))

delete_default_cube()

for area in bpy.context.screen.areas:
	if area.type == 'VIEW_3D':
		region = area.spaces[0].region_3d
		region.view_matrix = Matrix((
			( 0.4166,  0.9091, -0.0000,   -0.3976),
			(-0.5049,  0.2314,  0.8316,    3.8098),
			( 0.7560, -0.3464,  0.5554, -549.8699),
			( 0.0000,  0.0000,  0.0000,    1.0000)
			))

bpy.ops.wm.save_as_mainfile(filepath=blendpath)

)MULTILINE";
