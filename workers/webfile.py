import os
import yaml,json
from pathlib import Path
from celery import Celery
from celery.signals import worker_init, worker_process_init
import scanpy as sc
import pandas as pd
import json
import scipy
import csv
import os
from gzip import open as gzopen
import utils

app=Celery('webfile_task',broker='amqp://'+os.environ['RABBITMQ_HOST'],backend='redis://'+os.environ['REDIS_HOST'])

@worker_process_init.connect()
def on_worker_init(**_):
    print("WebFile Worker initiated")

@app.task(bind=True)
def task_list(self, *args, **kwargs):
    metafilename = Path(__file__).stem+".yml"
    taskobject = yaml.safe_load(open(metafilename,'r'))
    return taskobject

@app.task(bind=True)
def create_files(self, qcparams, **kwargs):
  self.update_state(state="STARTED")
  self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})
  config=utils.load_configuration()
  aws_s3=utils.AWS_S3(config)
  
  holder = []
  holder2 = []
  holder3 = []
  holder4 = []
  jsonStruct = {}
  jsonStruct2 = {}
  h5adPath = qcparams['aws_path']
  temp_dir = config['TEMP_DIRECTORY']
  path = Path(temp_dir).joinpath(h5adPath)
  RNA_flag = qcparams['rna_flag']

  if not RNA_flag:
    '''Begin to process the motif h5ad'''
    downloaded_filename_Motif = aws_s3.getFileObject('{}/motifs.h5ad'.format(h5adPath))
    adata=sc.read(downloaded_filename_Motif)
    if scipy.sparse.issparse(adata.X):
      adata.X = adata.X.toarray()
    df = pd.DataFrame(adata.X.transpose())
    f = gzopen('{}/motifSummation.txt.gz'.format(path), 'wt')
    f.write(str(adata.n_obs)+'\n')
    f.close()
    df.to_csv('{}/motifSummation.txt.gz'.format(path), float_format='%7.2f', index=False, header=False, sep=',', mode='a', compression='gzip')
    with gzopen('{}/motifNames.txt.gz'.format(path), 'wt') as employee_file2:
      for i in range(adata.n_vars):
        employee_file2.write(adata.var['mvp.variable'].index[i])
        employee_file2.write(',')
        
        
    adata.X = adata.X - adata.X.min() + 1
    adata.obs['clusters'] = adata.obs['clusters'].astype('category').values
    sc.tl.rank_genes_groups(adata, 'clusters', n_genes= 10, use_raw=False)
    holder = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    jsonStruct[1] = holder.values.tolist()
    
    adata2m=adata.copy()
    adata2m.X = -adata2m.X
    adata2m.X = adata2m.X - adata2m.X.min() + 1
    adata2m.obs['clusters'] = adata2m.obs['clusters'].astype('category').values
    sc.tl.rank_genes_groups(adata2m, 'clusters', n_genes= 10, use_raw=False)
    holder2 = pd.DataFrame(adata2m.uns['rank_genes_groups']['names'])
    jsonStruct[-1] = holder2.values.tolist()

    with open("{}/topTen_motifs.json".format(path), "w") as outfile:
        json.dump(jsonStruct, outfile)
        
  self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 50})
  
  '''Begin to process the gene h5ad'''
  downloaded_filename_Gene = aws_s3.getFileObject('{}/genes.h5ad'.format(h5adPath))
  adata2=sc.read(downloaded_filename_Gene)
  if scipy.sparse.issparse(adata2.X):
      adata2.X = adata2.X.toarray()
  df2 = pd.DataFrame(adata2.X.transpose())
  f2 = gzopen('{}/geneSummation.txt.gz'.format(path), 'wt')
  f2.write(str(adata2.n_obs)+'\n')
  f2.close()
  df2.to_csv('{}/geneSummation.txt.gz'.format(path), float_format='%7.2f', index=False, header=False, sep=',', mode='a', compression='gzip')

  with gzopen('{}/geneNames.txt.gz'.format(path), 'wt') as employee_file4:
    for i in range(adata2.n_vars):
      employee_file4.write(adata2.var['vst.variable'].index[i])
      employee_file4.write(',')
      
  with gzopen('{}/data.csv.gz'.format(path), 'wt') as employee_file5:
    write = csv.writer(employee_file5)
    for i in range(adata2.n_obs):
      if not RNA_flag:
        data = [adata2.obs['clusters'][i], adata2.obsm['spatial'][i].tolist(), adata2.obsm['X_umap'][i].tolist(), adata2.obs['TSSEnrichment'][i], adata2.obs['nFrags'][i]]
      else:
        data = [adata2.obs['clusters'][i], adata2.obsm['spatial'][i].tolist(), adata2.obsm['X_umap'][i].tolist(), adata2.obs['n_genes_by_counts'][i], adata2.obs['total_counts'][i]]
      write.writerow(data)


    adata2.X = adata2.X - adata2.X.min() + 1
    adata2.obs['clusters'] = adata2.obs['clusters'].astype('category').values
    sc.tl.rank_genes_groups(adata2, 'clusters', n_genes= 10, use_raw=False)
    holder3 = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    jsonStruct2[1] = holder3.values.tolist()
    
    adata2g=adata2.copy()
    adata2g.X = -adata2g.X
    adata2g.X = adata2g.X - adata2g.X.min() + 1
    adata2g.obs['clusters'] = adata2g.obs['clusters'].astype('category').values
    sc.tl.rank_genes_groups(adata2g, 'clusters', n_genes= 10, use_raw=False)
    holder4 = pd.DataFrame(adata2g.uns['rank_genes_groups']['names'])
    jsonStruct2[-1] = holder4.values.tolist()

    with open("{}/topTen_genes.json".format(path), "w") as outfile:
        json.dump(jsonStruct2, outfile)
    
    Path(downloaded_filename_Gene).unlink()
    Path(downloaded_filename_Motif).unlink()
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 100})
    
    

