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
import math
import pymysql
import time
import traceback
from natsort import natsorted

app=Celery('webfile_task',broker='amqp://'+os.environ['RABBITMQ_HOST'],backend='redis://'+os.environ['REDIS_HOST'])
@worker_process_init.connect()
def on_worker_init(**_):
    print("WebFile Worker initiated")

@app.task(bind=True)
def task_list(self, *args, **kwargs):
    metafilename = Path(__file__).stem+".yml"
    taskobject = yaml.safe_load(open(metafilename,'r'))
    return taskobject

def creatingOneThousand(totalNum):
  indexList = []
  numberOfSplits = divmod(totalNum, 1000)
  left = 0
  right = 0
  for i in range(numberOfSplits[0] + 1):
    if right >= totalNum: break
    if totalNum < 999: indexList.append((0, totalNum))
    else:
      left = i * 1000
      right = left + 999
      if right <= totalNum: indexList.append((left, right))
      else: indexList.append((left, totalNum))
  return indexList

def get_tissue_user_id(cursor, run_id, username):
  sql_get_userid = """ SELECT user_id FROM user_table WHERE username = %s"""
  sql_get_tissue_id = """ SELECT tissue_id from tissue_slides WHERE run_id = %s"""
  cursor.execute(sql_get_userid, (username,))
  user_id = cursor.fetchone()
  cursor.execute(sql_get_tissue_id, (run_id,))
  tissue_id = cursor.fetchone()
  if user_id:
    user_id = user_id[0]
  else:
    user_id = None
    raise Exception("User not found")
  if tissue_id:
    tissue_id = tissue_id[0]
  else:
    raise Exception("run id has no tissue id")
  return user_id, tissue_id

def write_join_table(cursor, job_id, tissue_id, db_connection):
  sql = """ INSERT INTO job_tissue_id_table (job_id, tissue_id) VALUES (%s, %s) """
  tup = (job_id, tissue_id)
  cursor.execute(sql, tup)
  db_connection.commit()

@app.task(bind=True)
def create_files(self, qcparams, **kwargs):
  db_connection = pymysql.connect(host=os.environ['MYSQL_HOST'], user=os.environ['MYSQL_USERNAME'], port = int(os.environ["MYSQL_PORT"]), password=os.environ['MYSQL_PASSWORD'], db=os.environ['MYSQL_DB'])
  username = kwargs.get('username', "")
  run_id = kwargs.get('run_id', "")
  description = kwargs.get("description", "")
  cursor = db_connection.cursor()
  user_id, tissue_id = get_tissue_user_id(cursor, run_id, username)
  current_time = int(time.time() * 1000)
  sql = """ INSERT INTO job_table (job_name, job_status, job_start_time, user_id, job_description) VALUES (%s, %s, %s, %s, %s) """
  tup = ("webfile.create_files", "INPROGRESS", current_time, user_id, description)
  cursor.execute(sql, tup)
  db_connection.commit()
  job_id = cursor.lastrowid
  write_join_table(cursor, job_id, tissue_id, db_connection)
  try:
    aws_bucket = qcparams.get('bucket_name', None)
    self.update_state(state="STARTED")
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})
    config=utils.load_configuration()
    print(config)
    if aws_bucket:
      config['AWS_BUCKET_NAME'] = aws_bucket
    aws_s3=utils.AWS_S3(config)
    
    holder2 = []
    holder3 = []
    holder4 = []
    jsonStruct = {}
    jsonStruct2 = {}
    h5adPath = qcparams['aws_path']
    temp_dir = config['TEMP_DIRECTORY']
    h5_index = h5adPath.rfind('/obj')
    path = Path(temp_dir).joinpath(h5adPath[:h5_index])
    RNA_flag = qcparams['rna_flag']
    
    # check that the path used in the next chunk exists
    if not os.path.exists('{}/summations'.format(path)):
      os.makedirs('{}/summations'.format(path))
    if not RNA_flag:
      '''Begin to process the motif h5ad'''
      downloaded_filename_Motif = aws_s3.getFileObject('{}/motifs.h5ad'.format(h5adPath))
      adata=sc.read(downloaded_filename_Motif)
      adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
      if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
      df = pd.DataFrame(adata.X.transpose())
      indexList = creatingOneThousand(adata.n_vars)
      for index in range(len(indexList)):
        sub = df.loc[indexList[index][0]:indexList[index][1]]
        f = gzopen('{}/summations/motifSummation{}.txt.gz'.format(path,index+1), 'wt')
        f.write(str(adata.n_obs)+'\n')
        f.close()
        sub.to_csv('{}/summations/motifSummation{}.txt.gz'.format(path,index+1), float_format='%6.2f', index=False, header=False, sep=',', mode='a', compression='gzip')
      with gzopen('{}/motifNames.txt.gz'.format(path), 'wt') as employee_file2:
        for i in range(adata.n_vars):
          employee_file2.write(adata.var['features'].index[i])
          employee_file2.write(',')
          
      adata.X = adata.X - adata.X.min() + 1
      adata.obs['clusters'] = adata.obs['clusters'].astype('category').values
      sc.tl.rank_genes_groups(adata, 'clusters', n_genes= 10, use_raw=False)
      holder = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
      nonSort  = list(holder.columns)
      clusters = natsorted(nonSort)
      jsonList = []
      for i in range(10):
        column = holder.loc[i]
        columnList = []
        for clust in natsorted(list(holder.columns)):
          value = column[clust]
          columnList.append(value)
        jsonList.append(columnList)
      jsonStruct[1] = jsonList
      
      adata2m=adata.copy()
      adata2m.X = -adata2m.X
      adata2m.X = adata2m.X - adata2m.X.min() + 1
      adata2m.obs['clusters'] = adata2m.obs['clusters'].astype('category').values
      sc.tl.rank_genes_groups(adata2m, 'clusters', n_genes= 10, use_raw=False)
      holder2 = pd.DataFrame(adata2m.uns['rank_genes_groups']['names'])
      jsonList = []
      for i in range(10):
        column = holder2.loc[i]
        columnList = []
        for clust in natsorted(list(holder2.columns)):
          value = column[clust]
          columnList.append(value)
        jsonList.append(columnList)
      jsonStruct[-1] = jsonList

      with open("{}/topTen_motifs.json".format(path), "w") as outfile:
          json.dump(jsonStruct, outfile)
          
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 50})
    
    '''Begin to process the gene h5ad'''
    downloaded_filename_Gene = aws_s3.getFileObject('{}/genes.h5ad'.format(h5adPath))
    adata2=sc.read(downloaded_filename_Gene)
    adata2.__dict__['_raw'].__dict__['_var'] = adata2.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    multiSample = 'Sample' in adata2.obs and 'Condition' in adata2.obs
    if scipy.sparse.issparse(adata2.X):
        adata2.X = adata2.X.toarray()
    df2 = pd.DataFrame(adata2.X.transpose())
    indexList = creatingOneThousand(adata2.n_vars)
    for index in range(len(indexList)):
      sub = df2.loc[indexList[index][0]:indexList[index][1]]
      f2 = gzopen('{}/summations/geneSummation{}.txt.gz'.format(path,index+1), 'wt')
      f2.write(str(adata2.n_obs)+'\n')
      f2.close()
      sub.to_csv('{}/summations/geneSummation{}.txt.gz'.format(path,index+1), float_format='%6.2f', index=False, header=False, sep=',', mode='a', compression='gzip')

    with gzopen('{}/geneNames.txt.gz'.format(path), 'wt') as employee_file4:
      for i in range(adata2.n_vars):
        employee_file4.write(adata2.var['features'].index[i])
        employee_file4.write(',')
        
    with gzopen('{}/data.csv.gz'.format(path), 'wt') as employee_file5:
      write = csv.writer(employee_file5)
      for i in range(adata2.n_obs):
        if not RNA_flag:
          if not multiSample:
            data = [adata2.obs['clusters'][i], adata2.obsm['spatial'][i][0].round(1), adata2.obsm['spatial'][i][1].round(1), adata2.obsm['X_umap'][i][0].round(2), adata2.obsm['X_umap'][i][1].round(2), adata2.obs['TSSEnrichment'][i].round(2), round(math.log10(adata2.obs['nFrags'][i]),2)]
          else:
            data = [adata2.obs['clusters'][i], adata2.obsm['spatial'][i][0].round(1), adata2.obsm['spatial'][i][1].round(1), adata2.obsm['X_umap'][i][0].round(2), adata2.obsm['X_umap'][i][1].round(2), adata2.obs['TSSEnrichment'][i].round(2), round(math.log10(adata2.obs['nFrags'][i]),2), adata2.obs['Sample'][i], adata2.obs['Condition'][i]]
        else:
          if not multiSample:
            if 'C' in str(adata2.obs['clusters'][i]):
              data = [adata2.obs['clusters'][i], adata2.obsm['spatial'][i][0].round(1), adata2.obsm['spatial'][i][1].round(1), adata2.obsm['X_umap'][i][0].round(2), adata2.obsm['X_umap'][i][1].round(2), adata2.obs['nFeature_Spatial'][i], adata2.obs['nCount_Spatial'][i]]
            else:
              data = ['C'+str(int(adata2.obs['clusters'][i]) + 1), adata2.obsm['spatial'][i][0].round(1), adata2.obsm['spatial'][i][1].round(1), adata2.obsm['X_umap'][i][0].round(2), adata2.obsm['X_umap'][i][1].round(2), adata2.obs['nFeature_Spatial'][i], adata2.obs['nCount_Spatial'][i]]
          else:
            if 'C' in str(adata2.obs['clusters'][i]):
              data = [adata2.obs['clusters'][i], adata2.obsm['spatial'][i][0].round(1), adata2.obsm['spatial'][i][1].round(1), adata2.obsm['X_umap'][i][0].round(2), adata2.obsm['X_umap'][i][1].round(2), adata2.obs['nFeature_Spatial'][i], adata2.obs['nCount_Spatial'][i], adata2.obs['Sample'][i], adata2.obs['Condition'][i]]
            else:
              data = ['C'+str(int(adata2.obs['clusters'][i]) + 1), adata2.obsm['spatial'][i][0].round(1), adata2.obsm['spatial'][i][1].round(1), adata2.obsm['X_umap'][i][0].round(2), adata2.obsm['X_umap'][i][1].round(2), adata2.obs['nFeature_Spatial'][i], adata2.obs['nCount_Spatial'][i], adata2.obs['Sample'][i], adata2.obs['Condition'][i]]
        write.writerow(data)

      adata2.X = adata2.X - adata2.X.min() + 1
      adata2.obs['clusters'] = adata2.obs['clusters'].astype('category').values
      sc.tl.rank_genes_groups(adata2, 'clusters', n_genes= 10, use_raw=False)
      holder3 = pd.DataFrame(adata2.uns['rank_genes_groups']['names'])
      jsonList = []
      for i in range(10):
        column = holder3.loc[i]
        columnList = []
        for clust in natsorted(list(holder3.columns)):
          value = column[clust]
          columnList.append(value)
        jsonList.append(columnList)
      jsonStruct2[1] = jsonList
      
      adata2g=adata2.copy()
      adata2g.X = -adata2g.X
      adata2g.X = adata2g.X - adata2g.X.min() + 1
      adata2g.obs['clusters'] = adata2g.obs['clusters'].astype('category').values
      sc.tl.rank_genes_groups(adata2g, 'clusters', n_genes= 10, use_raw=False)
      holder4 = pd.DataFrame(adata2g.uns['rank_genes_groups']['names'])
      jsonList = []
      for i in range(10):
        column = holder4.loc[i]
        columnList = []
        for clust in natsorted(list(holder4.columns)):
          value = column[clust]
          columnList.append(value)
        jsonList.append(columnList)
      jsonStruct2[-1] = jsonList

      with open("{}/topTen_genes.json".format(path), "w") as outfile:
          json.dump(jsonStruct2, outfile)

      Path(downloaded_filename_Gene).unlink()
      if not RNA_flag: Path(downloaded_filename_Motif).unlink()


    completion_time = int(time.time() * 1000)
    stop_inx = h5adPath.rfind('h5')
    db_path = "S3://{}/{}".format(aws_bucket, h5adPath[:stop_inx])
    sql_run_update = """UPDATE `results_metadata` SET `results_folder_path` = %s WHERE tissue_id = %s"""
    tup = (db_path, tissue_id)
    cursor.execute(sql_run_update, tup)
    db_connection.commit()

    sql_done = "UPDATE `job_table` SET `job_status` = 'SUCCESS', `job_completion_time` = %s, `job_execution_time` = %s WHERE `job_table`.`job_id` = %s"
    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 100})
  except Exception as e:
    print("There was an exception")
    #traceback the error and print it
    traceback.print_exc()
    print(e)
    completion_time = int(time.time() * 1000)
    sql_done = "UPDATE `job_table` SET `job_status` = 'FAILURE', `job_completion_time` = %s, `job_execution_time` = %s WHERE `job_table`.`job_id` = %s"
  finally:
    execution_time = completion_time - current_time
    tup = (completion_time, execution_time, job_id)
    cursor.execute(sql_done, tup)
    db_connection.commit()
    
    

