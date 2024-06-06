
import os
from threading import local
import traceback
import yaml,json,csv
import PIL
from PIL import Image
from pathlib import Path
import time
import re
import shutil
from anndata import AnnData
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
from matplotlib.image import imread
import json
from pathlib import Path
Image.MAX_IMAGE_PIXELS = None
from celery import Celery
from celery.signals import worker_init, worker_process_init

import numpy as np 
import utils
import cv2
import time
import math

app=Celery('atlasbrowser_task',broker='amqp://'+os.environ['RABBITMQ_HOST'],backend='redis://'+os.environ['REDIS_HOST'])

@worker_process_init.connect()
def on_worker_init(**_):
    print("AtlasBrowser Worker initiated")

@app.task(bind=True)
def task_list(self, *args, **kwargs):
    metafilename = Path(__file__).stem+".yml"
    taskobject = yaml.safe_load(open(metafilename,'r'))
    return taskobject

def rotate_image_no_cropping(img, degree):
    (h, w) = img.shape[:2]
    (cX, cY) = (w // 2, h // 2)
    # rotate our image by 45 degrees around the center of the image
    M = cv2.getRotationMatrix2D((cX, cY), degree, 1.0)
    abs_cos = abs(M[0,0]) 
    abs_sin = abs(M[0,1])
    bound_w = int(h * abs_sin + w * abs_cos)
    bound_h = int(h * abs_cos + w * abs_sin)
    M[0, 2] += bound_w/2 - cX
    M[1, 2] += bound_h/2 - cY
    rotated = cv2.warpAffine(img, M, (bound_w, bound_h))
    return rotated
@app.task(bind=True)
def generate_spatial(self, qcparams, **kwargs):
    self.update_state(state="STARTED")
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})
    config=utils.load_configuration()
    
    ############################################
    # Pull parameters from qcparams file
    
    temp_dir = config['TEMP_DIRECTORY'] 
    upload_list=[]
    root_dir_spatial = qcparams['root_dir_spatial']
    root_dir_bsa = qcparams["root_dir_bsa"]
    bucket_name_spatial = qcparams.get("bucket_name_spatial")
    bucket_name_bsa = qcparams.get("bucket_name_bsa")
    
    metadata = qcparams['metadata']
    oldFiles = qcparams['files']
    scalefactors = qcparams['scalefactors']
    run_id = qcparams['run_id']
    tixel_positions = qcparams['mask']
    crop_coordinates = qcparams['crop_area']
    orientation = qcparams['orientation']
    barcodes = qcparams.get('barcodes', 2)
    rotation = (int(orientation['rotation']) % 360)
    bsa_filename = qcparams['bsa_filename']
    barcodes = qcparams["barcode_list"]
    updating_existing = qcparams.get('updating_existing', False)
    
    
    config["S3_BUCKET_NAME"] =  bucket_name_spatial
    aws_s3_spatial = utils.AWS_S3(config)
    
    config_bsa = config.copy()
    config_bsa["S3_BUCKET_NAME"] = bucket_name_bsa
    aws_s3_bsa = utils.AWS_S3(config_bsa)

    # next_gen_barcodes = True
    
    # metadata["replaced_24_barcodes"] = next_gen_barcodes

    temp_path = Path(temp_dir).joinpath(root_dir_spatial, run_id)

    # Filtering the images that are in the current bsa directory to ensure they don't include contents of spatial folder
    allFiles = [i for i in oldFiles if ('.json' not in i and 'spatial' not in i)]
    
    ### output directories (S3)
    spatial_dir = Path(root_dir_spatial).joinpath(run_id, 'spatial')
    figure_dir = Path(root_dir_spatial).joinpath(run_id, 'spatial', 'figure')
    metadata_filename = spatial_dir.joinpath('metadata.json')
    scalefactors_filename = spatial_dir.joinpath('scalefactors_json.json')
    tissue_hires_image_filename = spatial_dir.joinpath('tissue_hires_image.png')
    tissue_lowres_image_filename = spatial_dir.joinpath('tissue_lowres_image.png')
    tissue_positions_filename = spatial_dir.joinpath('tissue_positions_list.csv')
    
    
    ### local temp directories
    local_spatial_dir = Path(temp_dir).joinpath(spatial_dir)
    local_figure_dir = Path(temp_dir).joinpath(figure_dir)
    # parents = True specifies that if not already existsing those parent directories are made.
    # exist_ok indicates to not re-make the folders if they already exist
    local_spatial_dir.mkdir(parents=True, exist_ok=True)
    local_figure_dir.mkdir(parents=True, exist_ok=True)

    
    ### read barcodes information 
    row_count = math.sqrt(len(barcodes))
    
    self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 20})
    ### save metadata & scalefactors
    local_metadata_filename = local_spatial_dir.joinpath('metadata.json')
    local_scalefactors_filename = local_spatial_dir.joinpath('scalefactors_json.json')
    json.dump(metadata, open(local_metadata_filename,'w'), indent=4,sort_keys=True)
    upload_list.append([local_metadata_filename,metadata_filename])
    # adding metadata and scalefactors to the list to be uploaded to S3 Bucket
    if not updating_existing:
        upload_list.append([local_scalefactors_filename,scalefactors_filename])
    ### load image from s3
        for i in allFiles:
            vals = i.split("/")
            name = vals[len(vals) - 1]
            
            if bsa_filename == i:
                bsa_path = Path(temp_dir).joinpath(i)
                bsa_original = Image.open(str(bsa_path))
                bsa_img_arr = np.array(bsa_original, np.uint8)
                if rotation != 0 :
                    bsa_img_arr = rotate_image_no_cropping(bsa_img_arr, rotation)
                
                try:
                    postB_img_arr = bsa_img_arr[:, :, 2]
                    postB_source = Image.fromarray(postB_img_arr)
                except:
                    postB_source = Image.fromarray(bsa_img_arr)
                
                bsa_source = Image.fromarray(bsa_img_arr)
                
            elif "flow" in i.lower() or "fix" in i.lower() or "dapi":
                path = str(figure_dir.joinpath(name))
                aws_s3_spatial.moveFile(bucket_name_bsa,bucket_name_spatial, i, path)

        self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 45})

        ## generate cropped images using crop parameters
        cropped_bsa = bsa_source.crop((crop_coordinates[0], crop_coordinates[1], crop_coordinates[2], crop_coordinates[3]))
        cropped_postB = postB_source.crop((crop_coordinates[0], crop_coordinates[1], crop_coordinates[2], crop_coordinates[3]))

        ## high resolution
        tempName_bsa = local_figure_dir.joinpath("postB_BSA.tif")
        tempName_postB = local_figure_dir.joinpath("postB.tif")
        s3Name_bsa = figure_dir.joinpath("postB_BSA.tif")
        s3Name_postB = figure_dir.joinpath("postB.tif")
        cropped_bsa.save(tempName_bsa.__str__())
        cropped_postB.save(tempName_postB.__str__())
        # adding figure folder images to upload list
        upload_list.append([tempName_postB, s3Name_postB])
        upload_list.append([tempName_bsa, s3Name_bsa])

        height = cropped_postB.height
        width = cropped_postB.width
        self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 65})
        if width > height:
            factorHigh = 2000/width
            factorLow = 600/width
            high_res = cropped_postB.resize((2000, int(height * factorHigh)), Image.ANTIALIAS)
            low_res = cropped_postB.resize((600, int(height * factorLow)), Image.ANTIALIAS)
        else:
            factorHigh = 2000/height
            factorLow = 600/height
            high_res = cropped_postB.resize((int(width*factorHigh), 2000), Image.ANTIALIAS)
            low_res = cropped_postB.resize((int(width*factorLow), 600), Image.ANTIALIAS)

        local_hires_image_path = local_spatial_dir.joinpath('tissue_hires_image.png')
        local_lowres_image_path = local_spatial_dir.joinpath('tissue_lowres_image.png')
        high_res.save(local_hires_image_path.__str__())
        low_res.save(local_lowres_image_path.__str__())
        scalefactors["tissue_hires_scalef"] = factorHigh
        scalefactors["tissue_lowres_scalef"] = factorLow
    
        json.dump(scalefactors, open(local_scalefactors_filename,'w'), indent=4,sort_keys=True)
        upload_list.append([local_hires_image_path, tissue_hires_image_filename])
        upload_list.append([local_lowres_image_path, tissue_lowres_image_filename])
        
    self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 75})
    ### generate tissue_positions_list.csv
    local_tissue_positions_filename= local_spatial_dir.joinpath('tissue_positions_list.csv')
    tissue_positions_list = []
    tixel_pos_list= [x['position'] for x  in tixel_positions]
    f=open(local_tissue_positions_filename, 'w')
    csvwriter = csv.writer(f, delimiter=',',escapechar=' ',quoting=csv.QUOTE_NONE)
    for idx, b in enumerate(barcodes):
        colidx = int(idx/row_count)
        rowidx = int(idx % row_count)
        keyindex = tixel_pos_list.index([rowidx,colidx])
        coord_x = int(round(tixel_positions[keyindex]['coordinates']['x']))
        coord_y = int(round(tixel_positions[keyindex]['coordinates']['y']))
        val = 0
        if tixel_positions[keyindex]['value'] : val = 1
        datarow = [b, val, rowidx, colidx, coord_y , coord_x ]
        tissue_positions_list.append(datarow)    
        csvwriter.writerow(datarow)
    f.close()
    self.update_state(state="PROGRESS", meta={"position": "Finishing" , "progress" : 80})
    upload_list.append([local_tissue_positions_filename, tissue_positions_filename])
    
    #Generating compressed spatial folder
    compressed_name = Path(temp_path).joinpath(f"{run_id}_compressed_spatial") #Create path to new file name, in the root/run_id dir
    compressed_output_key_spatial = Path(root_dir_spatial).joinpath(run_id, f"{run_id}_compressed_spatial.zip") #Creating the path to the uploaded to aws, including .zip extension
    shutil.make_archive( base_name=compressed_name, format='zip', root_dir = local_spatial_dir) #creating the zipped file
    upload_path_zipped_spatial = Path(str(compressed_name) + ".zip") #Added the .zip extension back to the path for aws upload
    upload_list.append([upload_path_zipped_spatial, compressed_output_key_spatial])
    
    for local_filename, output_key in upload_list:
        aws_s3_spatial.uploadFile(str(local_filename), str(output_key))
    
    self.update_state(state="PROGRESS", meta={"position": "Finished" , "progress" : 100})
    out=[list(map(str, x)) for x in upload_list]
    #print(len(out))
    return out

