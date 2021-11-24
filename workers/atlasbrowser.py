import os
import traceback
import yaml,json,csv
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
    run_id = qcparams['run_id']
    tixel_positions = qcparams['mask']
    crop_coordinates = qcparams['crop_area']
    orientation = qcparams['orientation']
    hflip = orientation['horizontal_flip']
    vflip = orientation['vertical_flip']
    rotation = int(orientation['rotation'])

    ### source image path
    image_dir = Path(root_dir).joinpath(run_id,'images')
    image_source_path = image_dir.joinpath('postB.tif')
    
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
    row_count = 50
    col_count = 50
    local_barcodes_filename = 'data/atlasbrowser/barcodes_50.txt'
    if len(tixel_positions) == 10000:
        local_barcodes_filename = 'data/atlasbrowser/barcodes_100.txt'
        row_count = 100
        col_count = 100

    barcodes = None
    with open(local_barcodes_filename,'r') as f:
        barcodes = f.read().splitlines()
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

    ### image flipping
    if hflip:
        source_image=cv2.flip(source_image, 1)
    if vflip:
        source_image=cv2.flip(source_image, 0)
    if rotation != 0 :
        for i in range(int(rotation/90)):
            source_image=cv2.rotate(source_image, cv2.ROTATE_90_CLOCKWISE)
    ### generate cropped images using crop parameters
    xs = [int(v['x']) for v in crop_coordinates]
    ys = [int(v['y']) for v in crop_coordinates]
    x1, y1 = [min(xs), min(ys)]
    x2, y2 = [max(xs), max(ys)]
        ## high resolution
    cropped_image = source_image[y1:y2,x1:x2]
    sf=scalefactors['tissue_hires_scalef']
    dim = [int(x*sf) for x in source_image.shape[0:2]]
    resized_image = cv2.resize(cropped_image, dim, interpolation=cv2.INTER_AREA)
    local_hires_image_path = local_spatial_dir.joinpath('tissue_hires_image.png')
    cv2.imwrite(local_hires_image_path.__str__(), resized_image, [cv2.IMWRITE_PNG_COMPRESSION, 9])
        ## low resolutions
    local_lowres_image_path = local_spatial_dir.joinpath('tissue_lowres_image.png')
    sf=scalefactors['tissue_lowres_scalef']
    dim = [int(x*sf) for x in source_image.shape[0:2]]
    resized_image = cv2.resize(cropped_image, dim, interpolation=cv2.INTER_AREA)
    cv2.imwrite(local_lowres_image_path.__str__(), resized_image, [cv2.IMWRITE_PNG_COMPRESSION, 9])
    upload_list.append([local_hires_image_path, tissue_hires_image_filename])
    upload_list.append([local_lowres_image_path, tissue_lowres_image_filename])

    ### generate tissue_positions_list.csv
    local_tissue_positions_filename= local_spatial_dir.joinpath('tissue_positions_list.csv')
    tissue_positions_list = []
    tixel_pos_list= [x['position'] for x  in tixel_positions]
    f=open(local_tissue_positions_filename, 'w')
    csvwriter = csv.writer(f, delimiter=',',escapechar=' ',quoting=csv.QUOTE_NONE)
    for idx, b in enumerate(barcodes):
        rowidx = int(idx/row_count)
        colidx = idx % row_count
        keyindex = tixel_pos_list.index([rowidx,colidx])
        coord_x = tixel_positions[keyindex]['coordinates']['x']
        coord_y = tixel_positions[keyindex]['coordinates']['y']
        val = 0
        if tixel_positions[keyindex]['value'] : val = 1
        datarow = [b, val, rowidx, colidx, coord_x , coord_y ]
        tissue_positions_list.append(datarow)    
        csvwriter.writerow(datarow)
    f.close()
    upload_list.append([local_tissue_positions_filename, tissue_positions_filename])

    ### concatenate tissue_positions_to gene expressions
    
    tissue_position_list_umi_genes_list = []
    for local_filename, output_key in upload_list:
        print("Copying {} to {}".format(local_filename, output_key))
        aws_s3.uploadFile(str(local_filename), str(output_key))

    self.update_state(state="PROGRESS", meta={"position": "Finishing" , "progress" : 100})
    out=[list(map(str, x)) for x in upload_list]
    print(len(out))
    return out


