from fileinput import filename
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

import utils

app=Celery('core_task',broker='amqp://'+os.environ['RABBITMQ_HOST'],backend='redis://'+os.environ['REDIS_HOST'])

@worker_process_init.connect()
def on_worker_init(**_):
    print("Creation Worker initiated")

@app.task(bind=True)
def task_list(self, *args, **kwargs):
    metafilename = Path(__file__).stem+".yml"
    taskobject = yaml.safe_load(open(metafilename,'r'))
    return taskobject

@app.task(bind=True)
def create_files(self, qcparams, **kwargs):
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)
    self.update_state(state="STARTED")
    ## config
    temp_dir = config['TEMP_DIRECTORY']
    upload_list=[]
    ## parameter parsing
    path = qcparams['path']
    data = qcparams['data']
    file_type = qcparams['file_type']
    file_name = qcparams['file_name']

    self.update_state(state="PROGRESS", meta={"position": "preparation" , "progress" : 0})

    path_dir = Path(path)
    data_filename = path_dir.joinpath(file_name)
    local_path_dir = Path(temp_dir).joinpath(path_dir)
    local_path_dir.mkdir(parents=True, exist_ok=True)
    local_data_filename = local_path_dir.joinpath(file_name)
    if file_type == 'json':
          json.dump(data, open(local_data_filename,'w'), indent=4,sort_keys=True)
          upload_list.append([local_data_filename,data_filename])
          self.update_state(state="PROGRESS", meta={"position": "running" , "progress" : 40})
          for local_filename, output_key in upload_list:
            aws_s3.uploadFile(str(local_filename), str(output_key))
            self.update_state(state="PROGRESS", meta={"position": "Finished" , "progress" : 100})
            out=[list(map(str, x)) for x in upload_list]
            return out
    

