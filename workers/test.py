import os
import traceback
import yaml,json
from pathlib import Path
import time

from celery import Celery
from celery.signals import worker_init, worker_process_init

import utils

app=Celery('core_task',broker='amqp://'+os.environ['RABBITMQ_HOST'],backend='redis://'+os.environ['REDIS_HOST'])

@worker_process_init.connect()
def on_worker_init(**_):
    print("Core Worker initiating")

@app.task(bind=True)
def task_list(self, *args, **kwargs):
    metafilename = Path(__file__).stem+".yml"
    taskobject = yaml.safe_load(open(metafilename,'r'))
    return taskobject

@app.task(bind=True)
def simple_task(self, *args, **kwargs):
    return " ".join(args)

@app.task(bind=True)
def long_task_fail(self, *args, **kwargs):
    self.update_state(state="STARTED")
    config=utils.load_configuration()
    limit=args[0]
    for i in range(limit):
        self.update_state(state="PROGRESS", meta={"position": "Process_{}".format(i) , "progress" : (i+1)/limit*100})
        time.sleep(1)
        if i>= limit*0.8: 
            raise Exception("Somehow failed")
    return "Successfully tested"

@app.task(bind=True)
def long_task_success(self, *args, **kwargs):
    self.update_state(state="STARTED")
    config=utils.load_configuration()
    limit=args[0]
    for i in range(limit):
        self.update_state(state="PROGRESS", meta={"position": "Process_{}".format(i) , "progress" : (i+1)/limit*100})
        time.sleep(1)
    return "Successfully tested"

@app.task(bind=True)
def rename_s3(self, *args, **kwargs):
    root_directory=args[0]
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)
    self.update_state(state="STARTED")
    filelist = aws_s3.getFileList(root_directory)
    filelist = list(map(lambda x: [x, x.replace('.out/','/out/')], filelist))
    for fn in filelist:
        src, dest = fn
        if src != dest:
            print("Copying {} to {}".format(src,dest))
            aws_s3.copyFile(src,dest)
            print("Deleting {}".format(src))
            aws_s3.deleteFile(src)
    return filelist;