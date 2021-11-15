import os
import traceback
import yaml,json
from pathlib import Path
import time

from celery import Celery
from celery.signals import worker_init, worker_process_init

app=Celery('core_task',broker='amqp://'+os.environ['RABBITMQ_HOST'],backend='redis://'+os.environ['REDIS_HOST'])

@worker_process_init.connect()
def on_worker_init(**_):
    print("Core Worker initiating")

@app.task(bind=True)
def test_task(self, *args, **kwargs):
    return " ".join(args)

@app.task(bind=True)
def test_long_task(self, *args, **kwargs):
    self.update_state(state="STARTED")
    limit=args[0]
    for i in range(limit):
        self.update_state(state="PROGRESS", meta={"position": "Process_{}".format(i) , "progress" : (i+1)/limit*100})
        time.sleep(1)
        if i>= limit*0.8: 
            raise Exception("Somehow failed")
    return "Successfully tested"
