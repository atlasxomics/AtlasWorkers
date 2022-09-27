import os
import traceback
import yaml,json
from pathlib import Path
import time

from celery import Celery
from celery.signals import worker_init, worker_process_init

import scanpy as sc
import numpy as np
import pandas as pd
import scipy

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
    filename,requested_genes, selected, rankGeneKey = args
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})
    downloaded_filename = aws_s3.getFileObject(filename)
    adata=sc.read(downloaded_filename)
    holder2 = []
    out={}
    if (len(selected) > 0):
        length = adata.n_obs
        lasso = ["unselected" for x in range(length)]
        for i in selected:
            lasso[i] = "selected"
        adata.obs["lasso"] = lasso
        if rankGeneKey == 1:
          sc.tl.rank_genes_groups(adata, 'lasso', n_genes= 10, use_raw=False)
          holder2 = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
        else:
          adata2=adata
          adata2.X = -adata2.X
          sc.tl.rank_genes_groups(adata, 'lasso', n_genes= 10, use_raw=False)
          holder2 = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
        out['top_selected'] = holder2['selected'].values.tolist()
        out['cluster_names'] = []
        out['top_ten'] = []
        out['clusters'] = []
        out['coordinates'] = []
        out['coordinates_umap'] = []
        out['genes'] = {}
        out['genes_summation'] = []
        out['total_genes'] = []
        out['TSS'] = []
        out['frags'] = []
        self.update_state(state="PROGRESS", meta={"position": "Finishing" , "progress" : 100})
        return out
    else:
        out['top_selected'] = holder2
        if rankGeneKey == 1:
          adata.obs['clusters'] = adata.obs['clusters'].astype('category').values
          sc.tl.rank_genes_groups(adata, 'clusters', n_genes= 10, use_raw=False)
          holder = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
        else:
          adata2=adata
          adata2.X = -adata2.X
          adata2.obs['clusters'] = adata2.obs['clusters'].astype('category').values
          sc.tl.rank_genes_groups(adata2, 'clusters', n_genes= 10, use_raw=False)
          holder = pd.DataFrame(adata2.uns['rank_genes_groups']['names'])
        out['cluster_names'] = list(holder.columns)
        out['top_ten'] = holder.values.tolist()
        out['clusters']=adata.obs['clusters'].tolist()
        out['coordinates']=adata.obsm['spatial'].tolist()
        out['coordinates_umap']=adata.obsm['X_umap'].tolist()
        out['genes']={}
        out['genes_summation']=np.zeros(len(out['coordinates']))
        out['total_genes'] = adata.var_names.to_list()
        out['TSS'] = adata.obs['TSSEnrichment'].tolist()
        out['frags'] = adata.obs['nFrags'].tolist()
        for g_exp in requested_genes:
            try:
                out['genes'][g_exp]= list(map(lambda x: x[0],adata[:,g_exp].X.todense()))
            except:
                out['genes'][g_exp]= list(map(lambda x: x[0],adata[:,g_exp].X))
        for k,v in out['genes'].items():
            out['genes_summation']+=np.array(v)
        out['genes_summation']=out['genes_summation'].tolist()
        self.update_state(state="PROGRESS", meta={"position": "Finishing" , "progress" : 100})
        return out

@app.task(bind=True)
def seq_logo(self, *args, **kwargs):
  config=utils.load_configuration()
  aws_s3=utils.AWS_S3(config)

  filename, id = args
  motif_csv = aws_s3.getFileObject(filename)
  position = id.index('-')
  motif_id = id[:position] + '_' + id[position+1:]

  motif_pwm = pd.read_csv(motif_csv)

  if motif_id not in motif_pwm['motif'].tolist():
      raise Exception("Motif not found")

  motif_pwm = motif_pwm[motif_pwm['motif'] == motif_id]
  motif_pwm = motif_pwm.dropna(axis = 1)

  bases = ['A', 'C', 'G', 'T']
  motif_pwm.insert(0, 'base', bases)

  positions = motif_pwm.columns.to_list()
  positions = [i for i in positions if i not in ['base', 'motif']]
  seqlogo_scores = [list(zip(motif_pwm['base'], motif_pwm[i])) for i in positions]


  return seqlogo_scores

@app.task(bind=True)
def compute_cell_type(self, *args, **kwargs):
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)
    self.update_state(state="STARTED")
    filename,marker_genes = args
    marker_genes = dict(marker_genes)
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})
    downloaded_filename = aws_s3.getFileObject(filename)
    adata=sc.read(downloaded_filename)
    numCluster = int(adata.obs['clusters'].value_counts().count())
    average = {}
    sort = {}
    cellTypeTable = {}
    allGenes = [item for sublist in marker_genes.values() for item in sublist]
    for i in allGenes:
      if i not in average.keys():
        average[i] = []
      for num in range(1, numCluster + 1):
        try:
          avg = float(np.mean(adata[adata.obs["clusters"] == 'C'+str(num), i].X))
          average[i].append(('C'+str(num), avg))
        except KeyError:
          pass
    for i,j in average.items():
      if len(j) > 0:
        holder = sorted(j, key = lambda x: x[1], reverse=True)
        sort[i] = holder
    if numCluster >= 2:
      for i,j in sort.items():
        for num in range(1):
          if int(scipy.stats.ttest_ind(adata[adata.obs["clusters"] == sort[i][num][0], i].X, adata[adata.obs["clusters"] == sort[i][num + 1][0], i].X, alternative="greater", trim=0.3)[0]) > 20:
            for cell in marker_genes.keys():
              if i in marker_genes[cell]:
                if i not in cellTypeTable.keys():
                  cellTypeTable[sort[i][num][0]] = []
                cellTypeTable[sort[i][num][0]].append(cell)
          else:
            for cell in marker_genes.keys():
              if i in marker_genes[cell]:
                if i not in cellTypeTable.keys():
                  cellTypeTable[sort[i][num][0]] = []
                cellTypeTable[sort[i][num][0]].append('Undefined')
    return cellTypeTable
