import os
import traceback
import yaml,json
from pathlib import Path
import time

from celery import Celery
from celery.signals import worker_init, worker_process_init

import numpy as np 
import utils
import cv2

app=Celery('atlasbrowser_task',broker='amqp://'+os.environ['RABBITMQ_HOST'],backend='redis://'+os.environ['REDIS_HOST'])

@worker_process_init.connect()
def on_worker_init(**_):
    print("AtlasBrowser Worker initiated")

@app.task(bind=True)
def generate_spatial(self, qcparams, **kwargs):
    self.update_state(state="STARTED")
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)
    ## config
    temp_dir = config['TEMP_DIRECTORY']
    upload_list=[]
    ## parameter parsing
    root_dir = qcparams['root_dir']
    metadata = qcparams['metadata']
    scalefactors = qcparams['scalefactors']
    run_id = metadata['run']
    tixel_positions = qcparams['tixel_positions']
    crop_coordinates = metadata['points']

    ### source image path
    image_dir = Path(root_dir).joinpath(run_id,'images')
    image_source_path = image_dir.joinpath('postB_BSA.tif')
    
    ### output directories (S3)
    spatial_dir = Path(root_dir).joinpath(run_id,'out','Gene','raw','spatial')
    raw_dir = Path(root_dir).joinpath(run_id,'out','Gene','raw')
    metadata_filename = spatial_dir.joinpath('metadata.json')
    scalefactors_filename = spatial_dir.joinpath('scalefactors_json.json')
    tissue_hires_image_filename = spatial_dir.joinpath('tissue_hires_image.png')
    tissue_lowres_image_filename = spatial_dir.joinpath('tissue_lowres_image.png')
    tissue_positions_filename = spatial_dir.joinpath('tissue_positions_list.csv')
    tissue_position_list_umi_filename = spatial_dir.joinpath('tissue_positions_list_log_UMI_Genes.csv')
    barcode_filename = raw_dir.joinpath('barcodes.tsv')
    ### local temp directories
    local_spatial_dir = Path(temp_dir).joinpath(spatial_dir)
    local_spatial_dir.mkdir(parents=True, exist_ok=True)

    ### read barcodes information 
    local_barcodes_filename = 'data/atlasbrowser/barcodes_50.txt'
    if len(tixel_positions) == 10000:
        local_barcodes_filename = 'data/atlasbrowser/barcodes_100.txt'

    barcodes = None
    with open(local_barcodes_filename,'r') as f:
        barcodes = f.readlines()
    upload_list.append([local_barcodes_filename, barcode_filename])

    ### save metadata & scalefactors
    local_metadata_filename = local_spatial_dir.joinpath('metadata.json')
    local_scalefactors_filename = local_spatial_dir.joinpath('scalefactors.json')
    json.dump(metadata, open(local_metadata_filename,'w'), indent=4,sort_keys=True)
    json.dump(scalefactors, open(local_scalefactors_filename,'w'), indent=4,sort_keys=True)
    upload_list.append([local_metadata_filename,metadata_filename])
    upload_list.append([local_scalefactors_filename,scalefactors_filename])

    ### load image from s3 
    local_image_path = aws_s3.getFileObject(str(image_source_path))
    source_image = cv2.imread(str(local_image_path), cv2.IMREAD_COLOR)

    ### generate cropped images using crop parameters
    xs = [int(v[1]) for v in list(enumerate(crop_coordinates)) if v[0] % 2 == 1]
    ys = [int(v[1]) for v in list(enumerate(crop_coordinates)) if v[0] % 2 == 0]
    x1, y1 = [min(xs), min(ys)]
    x2, y2 = [max(xs), max(ys)]
        ## high resolution
    cropped_image = source_image[y1:y2,x1:x2]
    sf=scalefactors['tissue_hires_scalef']
    dim = [int(x*sf) for x in cropped_image.shape[0:2]]
    resized_image = cv2.resize(cropped_image, dim, interpolation=cv2.INTER_AREA)
    local_hires_image_path = local_spatial_dir.joinpath('tissue_hires_image.png')
    cv2.imwrite(local_hires_image_path.__str__(), resized_image, [cv2.IMWRITE_PNG_COMPRESSION, 9])
        ## low resolutions
    local_lowres_image_path = local_spatial_dir.joinpath('tissue_lowres_image.png')
    sf=scalefactors['tissue_lowres_scalef']
    dim = [int(x*sf) for x in cropped_image.shape[0:2]]
    resized_image = cv2.resize(cropped_image, dim, interpolation=cv2.INTER_AREA)
    cv2.imwrite(local_lowres_image_path.__str__(), resized_image, [cv2.IMWRITE_PNG_COMPRESSION, 9])
    upload_list.append([local_hires_image_path, tissue_hires_image_filename])
    upload_list.append([local_lowres_image_path, tissue_lowres_image_filename])

    ### generate tissue_positions_list.csv
    
    tissue_positions_list = []
    
    ### concatenate tissue_positions_to gene expressions
    
    tissue_position_list_umi_genes_list = []

    self.update_state(state="PROGRESS", meta={"position": "Finishing" , "progress" : 100})
    out=[list(map(str, x)) for x in upload_list]
    print(len(out))
    return out


