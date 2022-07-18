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
    # hflip = orientation['horizontal_flip']
    # vflip = orientation['vertical_flip']
    rotation = int(orientation['rotation'])
    
    ### source image path
    allFiles = [i for i in oldFiles if '.json' not in i and 'spatial' not in i]
    
    ### output directories (S3)
    spatial_dir = Path(root_dir).joinpath(run_id,'images','spatial')
    figure_dir = Path(root_dir).joinpath(run_id,'images','spatial', 'figure')
    raw_dir = Path(root_dir).joinpath(run_id,'out','Gene','raw')
    metadata_filename = spatial_dir.joinpath('metadata.json')
    scalefactors_filename = spatial_dir.joinpath('scalefactors_json.json')
    tissue_hires_image_filename = spatial_dir.joinpath('tissue_hires_image.png')
    tissue_lowres_image_filename = spatial_dir.joinpath('tissue_lowres_image.png')
    tissue_positions_filename = spatial_dir.joinpath('tissue_positions_list.csv')
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

    self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 20})
    barcodes = None
    with open(local_barcodes_filename,'r') as f:
        barcodes = f.read().splitlines()
    ### save metadata & scalefactors
    local_metadata_filename = local_spatial_dir.joinpath('metadata.json')
    local_scalefactors_filename = local_spatial_dir.joinpath('scalefactors_json.json')
    json.dump(metadata, open(local_metadata_filename,'w'), indent=4,sort_keys=True)
    upload_list.append([local_metadata_filename,metadata_filename])
    upload_list.append([local_scalefactors_filename,scalefactors_filename])
    self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 40})
    ### load image from s3 
    for i in allFiles:
        if 'postb' in i.lower():
            local_image_path = aws_s3.getFileObject(str(i))
            source_image = Image.open(str(local_image_path))
            source_original = Image.open(str(local_image_path))
            save = source_original.save(str(local_image_path))
            temp = re.compile("(.+\/)(.+)")
            res = temp.search(i).groups() 
            fileName = res[1]
            ### image flipping
            # if hflip:
            #     horiz = source_image.transpose(PIL.Image.FLIP_TOP_BOTTOM)
            #     source_image = horiz
            # if vflip:
            #     vert = source_image.transpose(PIL.Image.FLIP_LEFT_RIGHT)
            #     source_image = vert
            if rotation != 0 :
                # if rotation == 90:
                #     rotation = 270
                #     rotate = source_image.rotate(rotation, expand=True)
                #     source_image = rotate
                # if rotation == 270:
                #     rotation = 90
                #     rotate = source_image.rotate(rotation, expand=True)
                #     source_image = rotate
                rotate = source_image.rotate(rotation, expand=True)
                source_image = rotate
            pillow_source_image = source_image
            ### generate cropped images using crop parameters
            cropped_image = pillow_source_image.crop((crop_coordinates[0], crop_coordinates[1], crop_coordinates[2], crop_coordinates[3]))
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

    spatial_dir = Path(root_dir).joinpath(run_id,'images','spatial')
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
    self.update_state(state="PROGRESS", meta={"position": "Finished" , "progress" : 100})
    out=[list(map(str, x)) for x in upload_list]
    return out
    