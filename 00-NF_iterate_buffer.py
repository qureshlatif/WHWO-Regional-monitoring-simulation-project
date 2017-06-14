## Generates 1km buffer applied to each forest ##
## All analyses and simulations were staged within these buffers ##

# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# NF_iterate_buffer.py
# Created on: 2016-01-26 11:57:51.00000
#   (generated by ArcGIS/ModelBuilder)
# Description: 
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy
import os
from arcpy import env

file_workspace = "T:\\FS\\RD\\RMRS\\Science\\WTE\\Research\\RMRS-WHWO\\Qs\\ArcGIS\\Occ_sims\\R6\\"
env.workspace = file_workspace

# Local variables
NFS = file_workspace + "Eastern_R6_forests_n83z10.shp"
NFS_Fieldname = "FORESTNUMB"
var_Buffer = "1000 Meters"

# Create NF folders if they don't already exist.
for i in range(0,len(NF_codes)):
    if not os.path.exists(file_workspace+"NF_"+NF_codes[i]):
        os.mkdir(file_workspace+"NF_"+NF_codes[i])

# Generate buffers
NF_codes = ["FROR","GPWA","OKWA","MAOR","OCOR","DEOR","WAOR","COWA","UMOR","MHOR"]
for i in range(0,len(NF_codes)):
    env.workspace = file_workspace+"NF_"+NF_codes[i]+"\\"
    arcpy.MakeFeatureLayer_management(NFS, "NF_lyr", '"FID" = ' + str(i))
    arcpy.Buffer_analysis("NF_lyr","NF_Buff1km1.shp",var_Buffer,"FULL","ROUND","ALL","","PLANAR")
    arcpy.Delete_management("NF_lyr")
