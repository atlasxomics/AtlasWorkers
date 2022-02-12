import os
import traceback
import yaml,json,csv
from PIL import Image
from pathlib import Path
import time
import re
import shutil

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
def task_list(self, *args, **kwargs):
    metafilename = Path(__file__).stem+".yml"
    taskobject = yaml.safe_load(open(metafilename,'r'))
    return taskobject

@app.task(bind=True)
def generate_spatial(self, qcparams, **kwargs):
    self.update_state(state="STARTED")
    self.update_state(state="PROGRESS", meta={"Prepartion": "one" , "progress" : 0})
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)
    ## config
    temp_dir = config['TEMP_DIRECTORY']
    upload_list=[]
    ## parameter parsing
    root_dir = qcparams['root_dir']
    metadata = qcparams['metadata']
    oldFiles = qcparams['files']
    scalefactors = qcparams['scalefactors']
    run_id = qcparams['run_id']
    tixel_positions = qcparams['mask']
    crop_coordinates = qcparams['crop_area']
    orientation = qcparams['orientation']
    barcodes = qcparams['barcodes']
    hflip = orientation['horizontal_flip']
    vflip = orientation['vertical_flip']
    rotation = int(orientation['rotation'])
    
    ### source image path
    allFiles = [i for i in oldFiles if '.json' not in i and 'spatial' not in i]
    
    ### output directories (S3)
    spatial_dir = Path(root_dir).joinpath(run_id, 'images', 'spatial')
    figure_dir = Path(root_dir).joinpath(run_id, 'images', 'spatial', 'figure')
    raw_dir = Path(root_dir).joinpath(run_id,'images','spatial')
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
    local_figure_dir = Path(temp_dir).joinpath(figure_dir)
    local_figure_dir.mkdir(parents=True, exist_ok=True)
    ### read barcodes information 
    row_count = 50
    local_barcodes_filename = 'data/atlasbrowser/barcodes_50.txt'
    if barcodes == 2:
        local_barcodes_filename = 'data/atlasbrowser/barcodesv2_50.txt'
    if len(tixel_positions) == 10000:
        local_barcodes_filename = 'data/atlasbrowser/barcodes_100.txt'
        row_count = 100

    self.update_state(state="PROGRESS", meta={"running": "two" , "progress" : 20})
    barcodes = None
    with open(local_barcodes_filename,'r') as f:
        barcodes = f.read().splitlines()
    upload_list.append([local_barcodes_filename, barcode_filename])
    ### save metadata & scalefactors
    local_metadata_filename = local_spatial_dir.joinpath('metadata.json')
    local_scalefactors_filename = local_spatial_dir.joinpath('scalefactors_json.json')
    json.dump(metadata, open(local_metadata_filename,'w'), indent=4,sort_keys=True)
    upload_list.append([local_metadata_filename,metadata_filename])
    upload_list.append([local_scalefactors_filename,scalefactors_filename])
    self.update_state(state="PROGRESS", meta={"position": "two" , "progress" : 40})
    ### load image from s3 
    for i in allFiles:
        if 'postb' in i.lower():
            local_image_path = aws_s3.getFileObject(str(i))
            source_image = cv2.imread(str(local_image_path), cv2.IMREAD_UNCHANGED)
            temp = re.compile("(.+\/)(.+)")
            res = temp.search(i).groups() 
            fileName = res[1]
            ### image flipping
            if hflip:
                source_image=cv2.flip(source_image, 1)
            if vflip:
                source_image=cv2.flip(source_image, 0)
            if rotation != 0 :
                for x in range(int(rotation/90)):
                    source_image=cv2.rotate(source_image, cv2.ROTATE_90_CLOCKWISE)
            cv2.imwrite(str(local_image_path), source_image)
            pillow_source_image = Image.open(str(local_image_path))
            ### generate cropped images using crop parameters
            cropped_image = pillow_source_image.crop((crop_coordinates[0]["x"], crop_coordinates[0]["y"], crop_coordinates[1]["x"], crop_coordinates[1]["y"]))
                ## high resolution
            tempName = local_figure_dir.joinpath(fileName)
            s3Name = figure_dir.joinpath(fileName)
            cropped_image.save(tempName.__str__())
            upload_list.append([tempName, s3Name])
            height = cropped_image.height
            width = cropped_image.width
            if 'bsa' not in i.lower():
                self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 65})
                if width > height:
                    factorHigh = 2000/width
                    factorLow = 600/width
                    high_res = cropped_image.resize((2000, int(height*factorHigh)), Image.ANTIALIAS)
                    low_res = cropped_image.resize((600, int(height*factorLow)), Image.ANTIALIAS)
                else:
                    factorHigh = 2000/height
                    factorLow = 600/height
                    high_res = cropped_image.resize((int(width*factorHigh), 2000), Image.ANTIALIAS)
                    low_res = cropped_image.resize((int(width*factorLow), 600), Image.ANTIALIAS)

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
    self.update_state(state="PROGRESS", meta={"position": "Finishing" , "progress" : 80})
    upload_list.append([local_tissue_positions_filename, tissue_positions_filename])

    ### concatenate tissue_positions_to gene expressions
    
    tissue_position_list_umi_genes_list = []
    for local_filename, output_key in upload_list:
        #print("Copying {} to {}".format(local_filename, output_key))
        aws_s3.uploadFile(str(local_filename), str(output_key))
    #shutil.make_archive('spatialcomp', 'zip', local_spatial_dir)
    self.update_state(state="PROGRESS", meta={"position": "Finished" , "progress" : 100})
    out=[list(map(str, x)) for x in upload_list]
    #print(len(out))
    return out
