import os
import traceback
import yaml,json
from pathlib import Path
import time

from celery import Celery
from celery.signals import worker_init, worker_process_init

import scanpy as sc
import numpy as np 

import utils

app=Celery('core_task',broker='amqp://'+os.environ['RABBITMQ_HOST'],backend='redis://'+os.environ['REDIS_HOST'])

@worker_process_init.connect()
def on_worker_init(**_):
    print("GeneMatrix Worker initiated")

@app.task(bind=True)
def task_list(self, *args, **kwargs):
    metafilename = Path(__file__).stem+".yml"
    taskobject = yaml.safe_load(open(metafilename,'r'))
    return taskobject

@app.task(bind=True)
def compute_qc(self, *args, **kwargs):
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)

    self.update_state(state="STARTED")
    filename,requested_genes = args
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})
    downloaded_filename = aws_s3.getFileObject(filename)
    computed_filename= Path(downloaded_filename).parent.joinpath(Path(downloaded_filename).stem+"_computed.h5ad").__str__()
    adata=None
    if Path(computed_filename).exists():
        adata=sc.read(computed_filename)
    else:      
        self.update_state(state="PROGRESS", meta={"position": "reading" , "progress" : 10})
        adata=sc.read(downloaded_filename)
        self.update_state(state="PROGRESS", meta={"position": "initial_qc" , "progress" : 20})
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        self.update_state(state="PROGRESS", meta={"position": "computing_pca" , "progress" : 30})
        sc.pp.pca(adata)
        self.update_state(state="PROGRESS", meta={"position": "clustering" , "progress" : 70})
        sc.pp.neighbors(adata)
        self.update_state(state="PROGRESS", meta={"position": "generating_umap" , "progress" : 80})
        sc.tl.umap(adata)
        self.update_state(state="PROGRESS", meta={"position": "further linalg" , "progress" : 90})
        sc.tl.leiden(adata,key_added="clusters")
        self.update_state(state="PROGRESS", meta={"position": "writing results" , "progress" : 95})
        adata.write(computed_filename)
    self.update_state(state="PROGRESS", meta={"position": "summarizing" , "progress" : 99})
    out={}
    out['clusters']=adata.obs['clusters'].tolist()
    out['coordinates']=adata.obsm['spatial'].tolist()
    out['coordinates_umap']=adata.obsm['X_umap'].tolist()
    out['genes']={}
    out['genes_summation']=np.zeros(len(out['coordinates']))
    for g_exp in requested_genes:
        try:
            out['genes'][g_exp]= list(map(lambda x: x[0],adata[:,g_exp].X.todense().tolist()))
        except:
            out['genes'][g_exp]= list(map(lambda x: x[0],adata[:,g_exp].X.tolist()))
    for k,v in out['genes'].items():
        out['genes_summation']+=np.array(v)
    out['genes_summation']=out['genes_summation'].tolist()
    self.update_state(state="PROGRESS", meta={"position": "Finishing" , "progress" : 100})
    return out

@app.task(bind=True)
def compute_qc_dev(self, *args, **kwargs):
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)

    self.update_state(state="STARTED")
    filename,requested_genes = args
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})
    downloaded_filename = aws_s3.getFileObject(filename)
    computed_filename= Path(downloaded_filename).parent.joinpath(Path(downloaded_filename).stem+"_computed.h5ad").__str__()
    adata=None
    if Path(computed_filename).exists():
        adata=sc.read(computed_filename)
    else:      
        self.update_state(state="PROGRESS", meta={"position": "reading" , "progress" : 10})
        adata=sc.read(downloaded_filename)
        self.update_state(state="PROGRESS", meta={"position": "initial_qc" , "progress" : 20})
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        self.update_state(state="PROGRESS", meta={"position": "computing_pca" , "progress" : 30})
        sc.pp.pca(adata)
        self.update_state(state="PROGRESS", meta={"position": "clustering" , "progress" : 70})
        sc.pp.neighbors(adata)
        self.update_state(state="PROGRESS", meta={"position": "generating_umap" , "progress" : 80})
        sc.tl.umap(adata)
        self.update_state(state="PROGRESS", meta={"position": "further linalg" , "progress" : 90})
        sc.tl.leiden(adata,key_added="clusters")
        self.update_state(state="PROGRESS", meta={"position": "writing results" , "progress" : 95})
        adata.write(computed_filename)
    self.update_state(state="PROGRESS", meta={"position": "summarizing" , "progress" : 99})
    out={}
    out['clusters']=list(map(lambda x: float(x)*0.5, adata.obs['clusters'].tolist()))
    out['coordinates']=adata.obsm['spatial'].tolist()
    out['coordinates_umap']=adata.obsm['X_umap'].tolist()
    out['genes']={}
    out['genes_summation']=np.zeros(len(out['coordinates']))
    for g_exp in requested_genes:
        try:
            out['genes'][g_exp]= list(map(lambda x: x[0],adata[:,g_exp].X.todense().tolist()))
        except:
            out['genes'][g_exp]= list(map(lambda x: x[0],adata[:,g_exp].X.tolist()))
    for k,v in out['genes'].items():
        out['genes_summation']+=np.array(v)
    out['genes_summation']=out['genes_summation'].tolist()
    self.update_state(state="PROGRESS", meta={"position": "Finishing" , "progress" : 100})
    return out

