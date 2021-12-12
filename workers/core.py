import os
import traceback
import yaml,json
from pathlib import Path
import time
from ftplib import FTP as FTP
from celery import Celery
from celery.signals import worker_init, worker_process_init
import threading
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
def list_novogene_batch(self, host, port, username, password , source_root='raw_data', output_root="data", **kwargs):
    self.update_state(state="STARTED")
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)
    tempDir = Path(config['TEMP_DIRECTORY'])
    res=None
    # print(host,port,username,password,filename,output_key)
    self.update_state(state="PROGRESS", meta={"position": "Connecting" , "progress" : 0})
    ftp=FTP()
    ftp.set_debuglevel(0)
    ftp.connect(host,port)
    ftp.login(username,password)
    res=ftp.nlst(source_root)
    res=list(filter(lambda x: '/D' in x,res))
    run_ids = list(map(lambda x: x.split('/')[1].split('_')[0], res))
    download_list=[]
    for idx, run_path in enumerate(res):
        run_files=ftp.nlst(run_path)
        run_files=list(filter(lambda x: 'fq.gz' in x, run_files))
        for fn in run_files:
            fn_idx = fn.split('.')[0].split('_')[-1]
            outname = "read{}.fq.gz".format(fn_idx)
            outkey = "{}/{}/sequences/{}".format(output_root,run_ids[idx],outname)
            obj = {
                'run_id' : run_ids[idx],
                'source_filename' : fn,
                'size' : ftp.size(fn),
                'output_key' : outkey,
                'output_exists' : aws_s3.checkFileExist(outkey)
            }
            download_list.append(obj)
    ftp.close()

    return download_list


@app.task(bind=True)
def download_from_ftp_to_s3(self, host, port, username, password , filename, output_key, **kwargs):
    self.update_state(state="STARTED")
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)
    tempDir = Path(config['TEMP_DIRECTORY'])
    res=None
    self.update_state(state="PROGRESS", meta={"position": "Connecting" , "progress" : 0})

    remotePath = filename
    tempFilepath=tempDir.joinpath('download_{}.tmp'.format(utils.get_uuid()))

    def downloadFile(host, port, username, password, fn, ofn):
        ftp = FTP(timeout=3600)
        ftp.connect(host,port)
        ftp.login(username,password)
        ftpReply = ftp.sendcmd('TYPE I');  ## binary mode
        print(ftpReply)
        size = ftp.size(fn)
        sock = ftp.transfercmd('RETR ' + fn)
        mb=1024*1024
        with open(ofn,'wb') as f:
            def background(f):
                while True:
                    block = sock.recv(1024*1024)
                    if not block:
                        break
                    f.write(block)
                sock.close()
            t = threading.Thread(target=background,args=[f])
            t.start()
            while t.is_alive():
                pos=-1
                try:
                    pos=f.tell()
                except:
                    pass
                pct = int(pos/size * 100)
                print("{} / {}".format(pos, size))
                self.update_state(state="PROGRESS", meta={"position": "{:.0f}/{:.0f}MB".format(pos/mb, size/mb) , "progress" : pct})
                time.sleep(1)
            ftp.close()
            pos=f.tell()
            print("{} / {}".format(pos, size))

    print("Transfer begins")
    downloadFile(host, port, username, password, remotePath, tempFilepath.__str__())
    print("Transfer completed")
    self.update_state(state="PROGRESS", meta={"position": "Copying" , "progress" : 100})
    print("File written")
    print("Copying {} to {}".format(tempFilepath.__str__(), output_key))
    aws_s3.uploadFile(tempFilepath.__str__(), str(output_key))
    print("Uploaded to S3")
    tempFilepath.unlink()
    print("Temporary file deleted")
    self.update_state(state="PROGRESS", meta={"position": "Finalizing" , "progress" : 100})
    time.sleep(3)
    return output_key

  

@app.task(bind=True)
def download_novogene_batch_to_s3(self, host, port, username, password , source_root='raw_data', output_root="data", **kwargs):
    self.update_state(state="STARTED")
    config=utils.load_configuration()
    aws_s3=utils.AWS_S3(config)
    tempDir = Path(config['TEMP_DIRECTORY'])
    res=None
    # print(host,port,username,password,filename,output_key)
    self.update_state(state="PROGRESS", meta={"position": "Connecting" , "progress" : 0})
    ftp=FTP()
    ftp.set_debuglevel(0)
    ftp.connect(host,port)
    ftp.login(username,password)
    res=ftp.nlst(source_root)
    res=list(filter(lambda x: '/D' in x,res))
    run_ids = list(map(lambda x: x.split('/')[1].split('_')[0], res))
    download_list=[]
    for idx, run_path in enumerate(res):
        run_files=ftp.nlst(run_path)
        run_files=list(filter(lambda x: 'fq.gz' in x, run_files))
        for fn in run_files:
            fn_idx = fn.split('.')[0].split('_')[-1]
            outname = "read{}.fq.gz".format(fn_idx)
            outkey = "{}/{}/sequences/{}".format(output_root,run_ids[idx],outname)
            obj = {
                'run_id' : run_ids[idx],
                'source_filename' : fn,
                'size': ftp.size(fn),
                'output_key' : outkey,
                'output_exists' : aws_s3.checkFileExist(outkey)
            }
            download_list.append(obj)
    ftp.close()

    global i
    global n_downloads
    n_downloads=len(download_list)
    for i, payload in enumerate(download_list):

        remotePath = payload['source_filename']
        tempFilepath=tempDir.joinpath('download_{}.tmp'.format(utils.get_uuid()))
        output_key = payload['output_key']
        def downloadFile(host, port, username, password, fn, ofn):
            global i
            global n_downloads
            ftp = FTP(timeout=3600)
            ftp.connect(host,port)
            ftp.login(username,password)
            ftpReply = ftp.sendcmd('TYPE I');  ## binary mode
            print(ftpReply)
            size = ftp.size(fn)
            sock = ftp.transfercmd('RETR ' + fn)
            mb=1024*1024
            with open(ofn,'wb') as f:
                def background(f):
                    while True:
                        block = sock.recv(1024*1024)
                        if not block:
                            break
                        f.write(block)
                    sock.close()
                t = threading.Thread(target=background,args=[f])
                t.start()
                while t.is_alive():
                    pos=-1
                    try:
                        pos=f.tell()
                    except:
                        pass
                    pct = int(pos/size * 100)
                    self.update_state(state="PROGRESS", meta={"position": "Downloading [{}/{}] {:.0f}mb/{:.0f}mb of {}, Dest : {}".format(i+1,n_downloads,pos/mb, size/mb, fn.split('/')[-1], str(output_key)) , "progress" : pct})
                    time.sleep(1)
                ftp.close()
                pos=f.tell()
                print("{} / {}".format(pos, size))

        print("Transfer begins")
        downloadFile(host, port, username, password, remotePath, tempFilepath.__str__())
        print("Transfer completed")
        self.update_state(state="PROGRESS", meta={"position": "Copying" , "progress" : 100})
        print("File written")
        print("Copying {} to {}".format(tempFilepath.__str__(), output_key))
        aws_s3.uploadFile(tempFilepath.__str__(), str(output_key))
        print("Uploaded to S3")
        tempFilepath.unlink()
        print("Temporary file deleted")
    self.update_state(state="PROGRESS", meta={"position": "Finalizing" , "progress" : 100})
    time.sleep(1)
    return res



