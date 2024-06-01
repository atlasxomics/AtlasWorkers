
import os
from threading import local
import yaml,json,csv
import PIL
from PIL import Image
from pathlib import Path
import shutil
import pandas as pd
from matplotlib.image import imread
import json
from pathlib import Path
Image.MAX_IMAGE_PIXELS = None
from celery import Celery
from celery.signals import worker_process_init
import utils
import cv2
import math
from modules.count_tixel import NucleiCounter

app=Celery('nuclei_task',broker='amqp://'+os.environ['RABBITMQ_HOST'],backend='redis://'+os.environ['REDIS_HOST'])

@worker_process_init.connect()
def on_worker_init(**_):
    print("Nuclei Worker initiated")

@app.task(bind=True)
def task_list(self, *args, **kwargs):
    metafilename = Path(__file__).stem+".yml"
    taskobject = yaml.safe_load(open(metafilename,'r'))
    return taskobject


@app.task(bind=True)
def generate_position_files(self, qcparams, **kwargs):
    self.update_state(state="STARTED")
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)
    
    ############################################
    # Pull parameters from qcparams file
    
    temp_dir = config['TEMP_DIRECTORY'] 
    upload_list=[]
    root_dir_spatial = qcparams['root_dir_spatial']
    bucket_name_spatial = qcparams.get("bucket_name_spatial")
    
    run_id = qcparams['run_id']
    tixel_positions = qcparams['mask']
    barcodes = qcparams["barcode_list"]
    dapi_path = qcparams["dapi_path"]
    tixel_width = qcparams["tixel_width"]
    user_given_count = qcparams["user_given"]
    thresh = qcparams['thresh']
    output_path = args['output_path_count']
    output_path_barcode = args['output_path_barcode']
    
    config["S3_BUCKET_NAME"] =  bucket_name_spatial
    

    temp_path = Path(temp_dir).joinpath(root_dir_spatial, run_id)

    ### output directories (S3)
    spatial_dir = Path(root_dir_spatial).joinpath(run_id, 'spatial')
    figure_dir = Path(root_dir_spatial).joinpath(run_id, 'spatial', 'figure')
    tissue_positions_filename = spatial_dir.joinpath('tissue_positions_list.csv')
    local_tissue_positions_filename = local_spatial_dir.joinpath('tissue_positions_list.csv')
    
    
    ### local temp directories
    local_spatial_dir = Path(temp_dir).joinpath(spatial_dir)
    local_figure_dir = Path(temp_dir).joinpath(figure_dir)
    local_spatial_dir.mkdir(parents=True, exist_ok=True)
    local_figure_dir.mkdir(parents=True, exist_ok=True)


    args = {'min_area_in_tixel':5, 'output_path':'/Users/joshuab/Desktop/Nuclei-Counting/go/googoo.csv'}
    image_path = dapi_path
    tissue_pos_path = aws_s3.getFileObject(tissue_positions_filename)
    channel = 2
    tixel_width_pixels = tixel_width
    
    min_area_in_tixel = args['min_area_in_tixel']
        
    
    column_names = ["barcode", "on_tissue", "row_inx", "col_inx", "y_coord", "x_coord"]
    dtypes = {
        "barcode": str,
        "on_tissue": str,
        "row_inx" : int,
        "col_inx" : int,
        "y_coord": int,
        "x_coord" : int,
    }
    
    tissue_postion_list = pd.read_csv(tissue_pos_path, sep=',', header=None, names=column_names, dtype=dtypes)
    
    image_object = aws_s3.getFileObject(dapi_path)
    read_image = cv2.imread(image_object.__str__(), cv2.IMREAD_COLOR)
    image = cv2.cvtColor(read_image, cv2.COLOR_BGR2GRAY)
    sel = int(thresh['c_value']) + 1
    sec = int(thresh['neighbor']) - 1
    if sel %2 == 0:
        sel+=1
    gray_image = cv2.adaptiveThreshold(image, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, sel, sec)
    
    self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 20})
    
    counting_obj = NucleiCounter(min_area_in_tixel=min_area_in_tixel)
    tissue_pos_cell_count = counting_obj.generate_tixel_mapping_cellpose(gray_image, tixel_width_pixels, tissue_postion_list)
    tissue_pos_cell_count.to_csv(output_path, index=False, header=False)
    
    qualify_barcodes = tissue_pos_cell_count[tissue_pos_cell_count['cell_count'] <= user_given_count]['barcode']
    qualify_barcodes.to_csv(output_path_barcode, index=False, header=False)

    
    '''
    OG CODE
    '''
    ### read barcodes information 
    row_count = math.sqrt(len(barcodes))
    
    self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 20})
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

