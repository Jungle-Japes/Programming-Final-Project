# -*- coding: utf-8 -*-
"""
Created on Tue May 07 18:51:42 2019

@author: Jo Belanger, Juan Cortes, Jacob Milley

"""

# import arcpy, check out spatial analyst license.
import arcpy, os
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")

# set environment workspace
# allow files to be overwritten
arcpy.env.workspace = r'S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\all_data'
arcpy.env.extent = r'S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\all_data\tinghir_dissolve.shp'
arcpy.env.overwriteOutput = True

#%%                 
            
# Rasterize four shapefile layers using the PolylineToRaster tool
# Output will be five rasters for each road type
road_list = ["primary_road.shp", 
		"secondary_road.shp", 
		"tertiary_road.shp",
		"unclassified_road.shp", 
		"footpath_road.shp"]

# Set output location
out_location = os.path.join(r'S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\all_data')

# Use a for loop to rasterize each of the road shapefiles
for road in road_list:   
    # use the FeatureToRaster conversion tool to rasterize each road.shp
    arcpy.FeatureToRaster_conversion(in_features= road,
                                      field = "type",
                                      out_raster =  os.path.join(out_location, '{}_raster.tif'.format(road)))

#%%
                              
# Reclassify each of the five rasters
# Output will be five reclassified rasters for each road  
# set workspace
arcpy.env.workspace = r"S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\all_data"
# list all rasters to the screen
print(arcpy.ListRasters())


#%% 

# reclassify primary roads
# use arcpy.Raster to make "primary" variable a raster object of that specific file.
primary = arcpy.Raster("primary_road.shp_raster.tif")

# reclassify primary roads
out_reclassify1 = Reclassify(in_raster = 'primary_road.shp_raster.tif',
                                # "Value" is the field that contains cell values
                                reclass_field = "Value",
                                # reclass all cell values of 1 to be 5
                                remap = RemapValue([[1, 1], [0, 100], ["NODATA", 100], [-9999, 100]]))                                                              
# set output location (note: different folder than your workspace)
# give your new file a name. (note: include .tif)
out_reclassify1.save(r'S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\all_data\primary_reclass.tif')

#%%

# use arcpy.Raster to make "secondary" variable a raster object of that specific file.
secondary = arcpy.Raster("secondary_road.shp_raster.tif")

# reclassify secondary roads
out_reclassify2 = Reclassify(in_raster = "secondary_road.shp_raster.tif",
                                # "Value" is the field that contains cell values
                                reclass_field = "Value",
                                # reclass all cell values of 1 to be 4
                                remap = RemapValue([[1, 2], [0, 100], ["NODATA", 100], [-9999, 100]]) )
# set output location
# name the file (note: include file extension)
out_reclassify2.save(r'S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\all_data\secondary_reclass.tif')

#%%

# reclassify tertiary roads
# use arcpy.Raster to make "tertiary" variable a raster object of that specific file.
tertiary = arcpy.Raster("tertiary_road.shp_raster.tif")

# reclassify tertiary roads
out_reclassify3 = Reclassify(in_raster = 'tertiary_road.shp_raster.tif',
                                # "Value" is the field that contains cell values
                                reclass_field = "Value",
                                # reclass all cell values of 1 to be 3
                                remap = RemapValue([[1, 3], [0, 100], ["NODATA", 100], [-9999, 100]]))
# set output location
# name the file (note: include file extension)
out_reclassify3.save(r'S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\all_data\tertiary_reclass.tif')

#%% 

# reclassify unclassified roads
# use arcpy.Raster to make "unclassified" variable a raster object of that specific file.
unclassified = "unclassified_road.shp_raster.tif"

# reclassify unclassified roads
out_reclassify4 = Reclassify(in_raster = "unclassified_road.shp_raster.tif",
                                # "Value" is the field that contains cell values
                                reclass_field = "Value",
                                # reclass all cell values of 1 to be 2
                                remap = RemapValue([[1, 4], [0, 100], ["NODATA", 100], [-9999, 100]]))
# set output location
# name the file (note: include file extension)
out_reclassify4.save(r'S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\all_data\unclassified_reclass.tif')

#%% 

# reclassify footpaths
# use arcpy.Raster to make "footpath" variable a raster object of that specific file.
footpath = "footpath_road.shp_raster.tif"
# reclassify footpaths to have a cell value of 1
out_reclassify5 = Reclassify(in_raster = footpath,
                                # "Value" is the field that contains cell values
                                reclass_field = "Value",
                                # reclass all cell values of 1 to be 1
                                remap = RemapValue([[1, 5], [0, 100], ["NODATA", 100], [-9999, 100]]))
# set output location
# name the file (note: include file extension)
out_reclassify5.save(r'S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\all_data\footpath_reclass.tif')


#%% 

# set workspace 
arcpy.env.workspace = r'S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\all_data'

raster1 = arcpy.Raster("primary_reclass.tif")
raster2 = arcpy.Raster("secondary_reclass.tif")
raster3 = arcpy.Raster("tertiary_reclass.tif")
raster4 = arcpy.Raster("unclassified_reclass.tif")
raster5 = arcpy.Raster("footpath_reclass.tif")
# create list of all 5 rasters
combine_raster = raster1 + raster2 + raster3 + raster4 + raster5
# join all five road rasters together
combine_raster.save("combine_raster.tif")
#raster_mosaic.save(r'S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\jo testing\output_data')

#%%

# begin Jacob & Juan's script

featureclass = "villages_clip.shp"
field_names = [f.name for f in arcpy.ListFields(featureclass)]

print field_names

#%%
            #COST DISTANCE RASTER!!!
count = 0
while count < 17:
    #arcpy.SearchCursor(in_table = "villages_clip.shp",
    #                   field_names = "FID")
    arcpy.MakeFeatureLayer_management(in_features = "ait_ouassif.shp",
                                      out_layer = "villages_clip_lyr",
                                      where_clause = "FID = " + str(count))
    outCostDistance = CostDistance(in_source_data = "villages_clip_lyr",
                                   in_cost_raster = 'combine_raster.tif')
    outCostDistance.save("out_cost_distance_" + str(count)+'.tif')               
    count = count + 1    
                   
#%%

village_raster_list = arcpy.ListRasters(wild_card = "out_cost_distance*", raster_type = "TIF")

for village_raster in village_raster_list:
    print(village_raster)


#%% #DO THE REST OF THESE STEPS 2 MORE TIMES (FOR SCHOOLS AND HEALTH CENTERS)

# Creating an empty list to store all of point shapefiles for each village after extracting value to points
village_community_scores_list = []
# This for loop iterates through the list of cost distance rasters
# It first fills the variable "out_village_community" with the rasters name and adds addtional information for the new name
# It then extracts values to point using the community centers shapefile as the in point features, the in raster from the list and the out point shapefile as specified above
# It then stores this new shapefile in the list created at the beginning

for village_raster in village_raster_list: 
    out_village_community = '{}.shp'.format(village_raster)
    new_out_village_community = out_village_community.split('.')[0]+'_comm.shp'
    ExtractValuesToPoints(in_point_features = "community_centers.shp",
                         in_raster = village_raster,
                         out_point_features = new_out_village_community)      
    village_community_scores_list.append(new_out_village_community)

village_community_scores_list.sort()

print village_community_scores_list

#%%
# The function file_sorter is supposed to take a shapefile, and sort it based on a specific field using the
# Sort_management tool. It takes two arguments, the original filename, and the new file name. 
# It then returns a brand new shapefile, with the sorted categorty.

def file_sorter(filename, new_filename):
    
    in_features = filename
    out_dataset = new_filename
    
    arcpy.Sort_management(in_features,
                      out_dataset,
                      sort_field = [["RASTERVALU", "ASCENDING"]])
   
    return out_dataset
    

#%%
# The min_finder function finds the minimum value of a field, and stores the first value (0 position) of a
# shapefile. That value is then returned is the minimum point of that shapefile. 

def min_finder(filename,field):
    
    in_features = filename
    arcpy.MakeFeatureLayer_management(in_features, out_layer = "mylayer")
   
# MakeFeatureLayer loads the file into memory, which allows it to be searched with a cursor.   
   
    field_names = [field]
     
    storage = []

# Establsihes the variables that will be used in later steps
     
    with arcpy.da.SearchCursor("mylayer", field_names) as cursor:
         for row in cursor:
             storage.append(row[0])

# the search cursor searches through the specified field within the specified layer, takes the first value within
# that field, and stores it with the storage container

    arcpy.Delete_management("mylayer")
    
    return storage[0]

# returns the stored value

#%%

# the min_creator function will make new shapefiles out of those previously sorted villages that contain only one point, the
# minimum value, with all of its attributes. It will require the creation of FieldDelimiters so that an SQL query can be run
# on it. CopyFeatures_management creates the new shapefile with the SQL applied to it.

def min_creator_comm(filename, field):
    
    in_features = filename
    field_names = field
    
    newField1 = arcpy.AddFieldDelimiters(in_features, field_names)
    
    storage = min_finder(filename,field)    
    
    min_SQL = newField1 + " = " + str(storage)
    
       # location = list1.index(file_a)    
    
       # new_min_point_layer = '{}_community_min'.format(filename)    
    
    arcpy.MakeFeatureLayer_management(in_features,
                                  out_layer = 'new_min_point_layer',
                                  where_clause = min_SQL )
    
    new_min_point_layer_shapefile = filename.split('.')[0]+'_community_min.shp'
    print(new_min_point_layer_shapefile)
    arcpy.CopyFeatures_management('new_min_point_layer', new_min_point_layer_shapefile)                                  

#%%
 
# the following for loop makes use of the file_sorter and the min_finder functions. First, all feature classes within
# a folder are listed,  then file sorter function sorts an attribute within those files in an ascending order. Then the
 # min_finder function takes the first object within a specified file, and stores it.
comm_file_list = arcpy.ListFeatureClasses("*comm.shp")

comm_file_list.sort()

for file_a in comm_file_list:
    print(file_a)
    file_sorter(filename = file_a, new_filename = "temp_sortedfile.shp")
    min_finder(filename = file_a ,field = "RASTERVALU")
    min_creator_comm(filename = file_a, field = "RASTERVALU")

final_comm_file_list = arcpy.ListFeatureClasses("*community_min.shp")

for file_a in final_comm_file_list:
    
    village_field = file_a.split('_')[3]
    
    arcpy.MakeFeatureLayer_management(in_features = file_a,
                                      out_layer = 'temp_file_field')   

    arcpy.AddField_management(in_table = 'temp_file_field',
                              field_name = "v_FID",
                              field_type = "FLOAT") 
                              
    with arcpy.da.UpdateCursor(in_table = 'temp_file_field',
                               field_names = "v_FID") as cursor:
                                   for row in cursor:
                                       arcpy.CalculateField_management(in_table = 'temp_file_field',
                                                                       field = "v_FID",
                                                                       expression = village_field)
#%% 
                                                                       
file_list = arcpy.ListFeatureClasses("*community_min.shp")

arcpy.Merge_management(inputs= file_list, output = "FINAL_community_min.shp")        

print "The script is completed"

#%% 

# Join Juan's merged layers to the village layer
# DO UP TO HERE 3 TIMES

# Creating an empty list to store all of point shapefiles for each village after extracting value to points
village_school_scores_list = []
# This for loop iterates through the list of cost distance rasters
# It first fills the variable "out_village_community" with the rasters name and adds addtional information for the new name
# It then extracts values to point using the community centers shapefile as the in point features, the in raster from the list and the out point shapefile as specified above
# It then stores this new shapefile in the list created at the beginning

for village_raster in village_raster_list: 
    out_village_school = '{}.shp'.format(village_raster)
    new_out_village_school = out_village_school.split('.')[0]+'_school.shp'
    ExtractValuesToPoints(in_point_features = "schools_clip.shp",
                         in_raster = village_raster,
                         out_point_features = new_out_village_school)      
    village_school_scores_list.append(new_out_village_school)

village_school_scores_list.sort()

print village_school_scores_list


#%%

# the min_creator function will make new shapefiles out of those previously sorted villages that contain only one point, the
# minimum value, with all of its attributes. It will require the creation of FieldDelimiters so that an SQL query can be run
# on it. CopyFeatures_management creates the new shapefile with the SQL applied to it.

def min_creator_school(filename, field):
    
    in_features = filename
    field_names = field
    
    newField1 = arcpy.AddFieldDelimiters(in_features, field_names)
    
    storage = min_finder(filename,field)    
    
    min_SQL = newField1 + " = " + str(storage)
    
       # location = list1.index(file_a)    
    
       # new_min_point_layer = '{}_community_min'.format(filename)    
    
    arcpy.MakeFeatureLayer_management(in_features,
                                  out_layer = 'new_min_point_layer',
                                  where_clause = min_SQL )
    
    new_min_point_layer_shapefile = filename.split('.')[0]+'_school_min.shp'
    print(new_min_point_layer_shapefile)
    arcpy.CopyFeatures_management('new_min_point_layer', new_min_point_layer_shapefile)       
#%%
 
# the following for loop makes use of the file_sorter and the min_finder functions. First, all feature classes within
# a folder are listed,  then file sorter function sorts an attribute within those files in an ascending order. Then the
 # min_finder function takes the first object within a specified file, and stores it.
school_file_list = arcpy.ListFeatureClasses("*school.shp")

school_file_list.sort()

for file_a in health_file_list:
    print(file_a)
    file_sorter(filename = file_a, new_filename = "temp_sortedfile.shp")
    min_finder(filename = file_a ,field = "RASTERVALU")
    min_creator_school(filename = file_a, field = "RASTERVALU")

final_school_file_list = arcpy.ListFeatureClasses("*school_min.shp")

for file_a in final_school_file_list:
    
    village_field = file_a.split('_')[3]
    
    arcpy.MakeFeatureLayer_management(in_features = file_a,
                                      out_layer = 'temp_file_field')   

    arcpy.AddField_management(in_table = 'temp_file_field',
                              field_name = "v_FID",
                              field_type = "FLOAT") 
                              
    with arcpy.da.UpdateCursor(in_table = 'temp_file_field',
                               field_names = "v_FID") as cursor:
                                   for row in cursor:
                                       arcpy.CalculateField_management(in_table = 'temp_file_field',
                                                                       field = "v_FID",
                                                                       expression = village_field)
    
file1_list = arcpy.ListFeatureClasses("*school_min.shp")

#%%
arcpy.Merge_management(inputs= file1_list, output = "FINAL_school_min.shp")        

print "The script is completed"

#%% #DO THE REST OF THESE STEPS 2 MORE TIMES (FOR SCHOOLS AND HEALTH CENTERS)

# Creating an empty list to store all of point shapefiles for each village after extracting value to points
village_health_scores_list = []
# This for loop iterates through the list of cost distance rasters
# It first fills the variable "out_village_community" with the rasters name and adds addtional information for the new name
# It then extracts values to point using the community centers shapefile as the in point features, the in raster from the list and the out point shapefile as specified above
# It then stores this new shapefile in the list created at the beginning

for village_raster in village_raster_list: 
    out_village_health = '{}.shp'.format(village_raster)
    new_out_village_health = out_village_health.split('.')[0]+'_health.shp'
    ExtractValuesToPoints(in_point_features = "health_centers.shp",
                         in_raster = village_raster,
                         out_point_features = new_out_village_health)      
    village_health_scores_list.append(new_out_village_health)

village_health_scores_list.sort()

print village_health_scores_list


#%%

# the min_creator function will make new shapefiles out of those previously sorted villages that contain only one point, the
# minimum value, with all of its attributes. It will require the creation of FieldDelimiters so that an SQL query can be run
# on it. CopyFeatures_management creates the new shapefile with the SQL applied to it.

def min_creator_health(filename, field):
    
    in_features = filename
    field_names = field
    
    newField1 = arcpy.AddFieldDelimiters(in_features, field_names)
    
    storage = min_finder(filename,field)    
    
    min_SQL = newField1 + " = " + str(storage)
    
       # location = list1.index(file_a)    
    
       # new_min_point_layer = '{}_community_min'.format(filename)    
    
    arcpy.MakeFeatureLayer_management(in_features,
                                  out_layer = 'new_min_point_layer',
                                  where_clause = min_SQL )
    
    new_min_point_layer_shapefile = filename.split('.')[0]+'_health_min.shp'
    print(new_min_point_layer_shapefile)
    arcpy.CopyFeatures_management('new_min_point_layer', new_min_point_layer_shapefile)       
#%%
 
# the following for loop makes use of the file_sorter and the min_finder functions. First, all feature classes within
# a folder are listed,  then file sorter function sorts an attribute within those files in an ascending order. Then the
 # min_finder function takes the first object within a specified file, and stores it.
health_file_list = arcpy.ListFeatureClasses("*health.shp")

health_file_list.sort()

for file_a in health_file_list:
    print(file_a)
    file_sorter(filename = file_a, new_filename = "temp_sortedfile.shp")
    min_finder(filename = file_a ,field = "RASTERVALU")
    min_creator_health(filename = file_a, field = "RASTERVALU")

final_health_file_list = arcpy.ListFeatureClasses("*health_min.shp")

for file_a in final_health_file_list:
    
    village_field = file_a.split('_')[3]
    
    arcpy.MakeFeatureLayer_management(in_features = file_a,
                                      out_layer = 'temp_file_field')   

    arcpy.AddField_management(in_table = 'temp_file_field',
                              field_name = "v_FID",
                              field_type = "FLOAT") 
                              
    with arcpy.da.UpdateCursor(in_table = 'temp_file_field',
                               field_names = "v_FID") as cursor:
                                   for row in cursor:
                                       arcpy.CalculateField_management(in_table = 'temp_file_field',
                                                                       field = "v_FID",
                                                                       expression = village_field)
    
file2_list = arcpy.ListFeatureClasses("*health_min.shp")

#%%
arcpy.Merge_management(inputs= file2_list, output = "FINAL_health_min.shp")        

print "The script is completed"
#%%

#JOIN x 3

arcpy.JoinField_management(in_data = "ait_ouassif.shp",
                         in_field = "FID",
                         join_table = "FINAL_community_min.shp",
                         join_field = "v_FID",
                         fields = ["LOCATION", "ASSOC_NAME", "ACTIVITY_T", "v_FID", "RASTERVALU"])
#%%                         
arcpy.JoinField_management(in_data = "ait_ouassif.shp",
                         in_field = "FID",
                         join_table = "FINAL_school_min.shp",
                         join_field = "v_FID",
                         fields = ["NOM_ecole", "Cycle", "v_FID", "RASTERVALU"])
#%%                         
arcpy.JoinField_management(in_data = "ait_ouassif.shp",
                         in_field = "FID",
                         join_table = "FINAL_health_min.shp",
                         join_field = "v_FID",
                         fields = ["NOM_ETABL", "v_FID", "RASTERVALU"]) 
                        
#%%

arcpy.AddField_management(in_table = "ait_ouassif.shp",
                          field_name = "School_W", 
                          field_type = "FLOAT")
                          
# Adding a new field for the weighted Health center score
arcpy.AddField_management(in_table = "ait_ouassif.shp", 
                          field_name = "Health_W", 
                          field_type = "FLOAT")
 
# Adding a new field for the weighted community center score
arcpy.AddField_management(in_table = "ait_ouassif.shp", 
                          field_name = "Comm_W", 
                          field_type = "FLOAT")

#%%

# Calculating a the new field using a predetermined weight and multiplying this times the original school score field
arcpy.CalculateField_management(in_table = "ait_ouassif.shp",
                                field = "School_W", 
                                expression = "!RASTERVA_2! * 0.3",
                                expression_type = "PYTHON")

# Calculating a the new field using a predetermined weight and multiplying this times the original health center score field
arcpy.CalculateField_management(in_table = "ait_ouassif.shp",
                                field = "Health_W", 
                                expression = "!RASTERVA_3! * 0.5",
                                expression_type = "PYTHON")

# Calculating a the new field using a predetermined weight and multiplying this times the original community center score field
arcpy.CalculateField_management(in_table = "ait_ouassif.shp",
                                field = "Comm_W", 
                                expression = "!RASTERVA_1! * 0.2",
                                expression_type = "PYTHON")

#%% 
 
# Adding a new field for the total accessability scores
arcpy.AddField_management(in_table = "ait_ouassif.shp",
                          field_name = "A_score", 
                          field_type = "FLOAT")
 
# Calculating this new field by adding the weighted scores for schools, health centers, and community centers
arcpy.CalculateField_management(in_table = "ait_ouassif.shp",
                                field = "A_score", 
                                expression = "!School_W! + !Health_W! + !Comm_W!",
                                expression_type = "PYTHON")
                                
#%% 
# Sorting this new field in ascending score because the lowest score is has the greatest accessability
out_dataset = r'S:\GEOG_4308.80_ProgrammingforGeospatialApp_201901\Class-Shared\Group_4\all_data'

arcpy.Sort_management(in_dataset = "ait_ouassif.shp", 
                      out_dataset = "A_score_sort",
                      sort_field = [["A_score", "DESCENDING"]])
#%%

in_table = "A_score_sort.shp"

field_names = ["DOUAR", "A_score"]

with arcpy.da.SearchCursor(in_table, field_names) as cursor:
    for row in cursor:
        print ('{0} : {1}'.format(row[0], row[1]))

#print final_list

print "The script is completed!!!"
