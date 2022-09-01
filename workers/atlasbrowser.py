
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
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})
    config=utils.load_configuration()
    bucket_name = config["S3_BUCKET_NAME"]
    if "bucket" in qcparams.keys():
        bucket_name = qcparams["bucket"]
        config["S3_BUCKET_NAME"] = bucket_name
    print(config)
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
    rotation = int(orientation['rotation'])

    #remove all files from the temp folder. To allevaite bugs being caused by figure folder being generated using old images.
    temp_path = Path(temp_dir).joinpath(root_dir, run_id)
    if os.path.exists(temp_path):
        for filename in os.listdir(temp_path):
            file_path = os.path.join(temp_path, filename)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete ' + file_path + " due to " + e)

    ### source image path
    allFiles = [i for i in oldFiles if '.json' not in i and 'spatial' not in i]
    ### output directories (S3)
    spatial_dir = Path(root_dir).joinpath(run_id, 'spatial')
    figure_dir = Path(root_dir).joinpath(run_id, 'spatial', 'figure')
    raw_dir = Path(root_dir).joinpath(run_id,'out','Gene','raw')
    metadata_filename = spatial_dir.joinpath('metadata.json')
    scalefactors_filename = spatial_dir.joinpath('scalefactors_json.json')
    tissue_hires_image_filename = spatial_dir.joinpath('tissue_hires_image.png')
    tissue_lowres_image_filename = spatial_dir.joinpath('tissue_lowres_image.png')
    tissue_positions_filename = spatial_dir.joinpath('tissue_positions_list.csv')
    ### local temp directories
    local_spatial_dir = Path(temp_dir).joinpath(spatial_dir)
    # parents = True specifies that if not already existsing those parent directories are made.
    # exist_ok indicates to not re-make the folders if they already exist
    local_spatial_dir.mkdir(parents=True, exist_ok=True)
    local_figure_dir = Path(temp_dir).joinpath(figure_dir)
    local_figure_dir.mkdir(parents=True, exist_ok=True)
    
    ### read barcodes information 
    row_count = 50
    local_barcodes_filename = 'data/atlasbrowser/bc50v1.txt'
    if barcodes == 2:
        local_barcodes_filename = 'data/atlasbrowser/bc50v2.txt'
    elif barcodes == 3:
        local_barcodes_filename = 'data/atlasbrowser/bc50v3.txt'
    elif barcodes == 4:
        local_barcodes_filename = 'data/atlasbrowser/bc50v4.txt'

    self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 20})
    barcodes = None
    with open(local_barcodes_filename,'r') as f:
        barcodes = f.read().splitlines()
    ### save metadata & scalefactors
    local_metadata_filename = local_spatial_dir.joinpath('metadata.json')
    local_scalefactors_filename = local_spatial_dir.joinpath('scalefactors_json.json')
    json.dump(metadata, open(local_metadata_filename,'w'), indent=4,sort_keys=True)
    # adding metadata and scalefactors to the list to be uploaded to S3 Bucket
    upload_list.append([local_metadata_filename,metadata_filename])
    upload_list.append([local_scalefactors_filename,scalefactors_filename])
    ### load image from s3
    for i in allFiles:
        vals = i.split("/")
        name = vals[len(vals) - 1]
        if "flow" in i.lower() or "fix" in i.lower():
            path = str(figure_dir.joinpath(name))
            print("old: " + i)
            print("new: " + path)
            aws_s3.moveFile(bucket_name, i, path)
        if 'postb_bsa' in i.lower():
            local_image_path = aws_s3.getFileObject(str(i))
            bsa_original = Image.open(str(local_image_path))
            bsa_source = bsa_original
            bsa_original.save(str(local_image_path))
            img_arr = np.array(bsa_original, np.uint8)
            
            postB_img_arr = img_arr[:, :, 2]
            postB_original = Image.fromarray(postB_img_arr)
            postB_source = postB_original
            postB_original.save(str(local_image_path))

            self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 45})

    
    if rotation != 0 :
        rotate_bsa = bsa_source.rotate(rotation, expand = False)
        bsa_source = rotate_bsa
        rotate_postB = postB_source.rotate(rotation, expand=False)
        postB_source = rotate_postB

    ### generate cropped images using crop parameters
    cropped_bsa = bsa_source.crop((crop_coordinates[0], crop_coordinates[1], crop_coordinates[2], crop_coordinates[3]))
    cropped_postB = postB_source.crop((crop_coordinates[0], crop_coordinates[1], crop_coordinates[2], crop_coordinates[3]))
    ## high resolution
    tempName_bsa = local_figure_dir.joinpath("postB_BSA.tif")
    tempName_postB = local_figure_dir.joinpath("postB.tif")
    s3Name_bsa = figure_dir.joinpath("postB_BSA.tif")
    s3Name_postB = figure_dir.joinpath("postB.tif")
    print("saving cropped image")
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
        # if os.path.exists(str(local_filename)):
        #   os.remove(str(local_filename))
    #shutil.make_archive('spatialcomp', 'zip', local_spatial_dir)
    self.update_state(state="PROGRESS", meta={"position": "Finished" , "progress" : 100})
    out=[list(map(str, x)) for x in upload_list]
    #print(len(out))
    return out


@app.task(bind=True)
def generate_h5ad(self, qcparams, **kwargs):
    self.update_state(state="STARTED")
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)
    #config
    temp_dir = config['TEMP_DIRECTORY']
    upload_list = []
    #parameters
    root_dir = qcparams['root_dir']
    run_id = qcparams['run_id']

    bar = Path(root_dir + '/'+ run_id + '/out/Gene/raw/barcodes.tsv')
    pathBar = aws_s3.getFileObject(str(bar))
    gene = Path(root_dir + '/'+ run_id + '/out/Gene/raw/features.tsv')
    pathGene = aws_s3.getFileObject(str(gene))
    matrix = Path(root_dir + '/'+ run_id + '/out/Gene/raw/matrix.mtx')
    pathMatrix = aws_s3.getFileObject(str(matrix))

    spatial_dir = Path(root_dir).joinpath(run_id,'spatial')
    raw_dir = Path(root_dir).joinpath(run_id,'out','Gene','raw')
    h5_dir = Path(root_dir).joinpath(run_id,'h5','obj')
    barcode_filename = raw_dir.joinpath('barcodes.tsv')
    gene_filename = raw_dir.joinpath('features.tsv')
    matrix_filename = raw_dir.joinpath('matrix.mtx')
    h5_filename = h5_dir.joinpath(f'{run_id}.h5ad')
    ### local temp directories
    local_spatial_dir = Path(temp_dir).joinpath(spatial_dir)
    local_raw_dir = Path(temp_dir).joinpath(raw_dir)
    local_raw_dir.mkdir(parents=True, exist_ok=True)
    local_h5_dir = Path(temp_dir).joinpath(h5_dir)
    local_h5_dir.mkdir(parents=True, exist_ok=True)

    local_barcode_filename = local_raw_dir.joinpath('barcodes.tsv')
    local_gene_filename = local_raw_dir.joinpath('features.tsv')
    local_matrix_filename = local_raw_dir.joinpath('matrix.mtx')
    local_h5ad_filename = local_h5_dir.joinpath(f'{run_id}.h5ad')
    upload_list.append([local_h5ad_filename, h5_filename])
    upload_list.append([local_barcode_filename, barcode_filename])
    upload_list.append([local_gene_filename, gene_filename])
    upload_list.append([local_matrix_filename, matrix_filename])

    f = open(pathBar, 'r')
    f.close()
    g = open(pathGene, 'r')
    g.close()
    m = open(pathMatrix, 'r')
    m.close()

    path = Path(root_dir + '/'+ run_id + '/out/Gene/raw')
    path2 = aws_s3.getFileObject(str(path))
    adata = sc.read_10x_mtx(str(path2), make_unique="true", var_names="gene_symbols")
    adata.uns["spatial"] = dict()
    adata.uns["spatial"]["0"] = dict()
    
    files = dict(
        tissue_positions_file = local_spatial_dir / 'tissue_positions_list.csv',
        scalefactors_json_file = local_spatial_dir / 'scalefactors_json.json',
        hires_image = local_spatial_dir / 'tissue_hires_image.png',
        lowres_image = local_spatial_dir / 'tissue_lowres_image.png'
    )
    adata.uns["spatial"]["0"]['images'] = dict()
    for res in ['hires', 'lowres']:
        try:
            adata.uns["spatial"]["0"]['images'][res] = imread(
                str(files[f'{res}_image'])
                )
        except Exception:
            raise OSError(f"Could not find {res}_image'")
    self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 40})
    # read json scalefactors
    print('good ole chum')
    adata.uns["spatial"]["0"]['scalefactors'] = json.loads(
        files['scalefactors_json_file'].read_bytes()
    )
    # read coordinates
    positions = pd.read_csv(files['tissue_positions_file'], header=None)
    positions.columns = [
        'barcode',
        'in_tissue',
        'array_row',
        'array_col',
        'pxl_col_in_fullres',
        'pxl_row_in_fullres',
    ]
        
    positions.index = positions['barcode']

    adata.obs = adata.obs.join(positions, how="left")

    adata.obsm['spatial'] = adata.obs[
        ['pxl_row_in_fullres', 'pxl_col_in_fullres']
    ].to_numpy()
            
    adata.obs.drop(
        columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'],
        inplace=True,
    )
    adata=adata[adata.obs['in_tissue']==1]
    self.update_state(state="PROGRESS", meta={"position": "Finishing" , "progress" : 80})
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith(('mt-','MT-', 'Mt-'))  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 30, :]
    sc.pp.normalize_total(adata, target_sum= 1e6, exclude_highly_expressed=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25)
    sc.tl.umap(adata, n_components=3)
    sc.tl.leiden(adata, key_added="clusters")
    adata.write(local_h5ad_filename)
    for local_filename, output_key in upload_list:
        aws_s3.uploadFile(str(local_filename), str(output_key))
        if os.path.exists(str(local_filename)):
            os.remove(str(local_filename))
    self.update_state(state="PROGRESS", meta={"position": "Finished" , "progress" : 100})
    out=[list(map(str, x)) for x in upload_list]
    return out
    