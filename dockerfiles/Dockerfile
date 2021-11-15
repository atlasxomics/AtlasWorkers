FROM ubuntu:20.04
WORKDIR /
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y &&\
    apt-get install build-essential python3-pip python3-venv libgtk2.0-dev -y &&\
    pip3 install uwsgi flask flask_restful flask-cors pymongo requests boto3 requests celery celery[redis] pika redis &&\
    pip3 install flask_jwt_extended uwsgi &&\
    pip3 install pyyaml warrant numpy==1.20 scanpy leidenalg


