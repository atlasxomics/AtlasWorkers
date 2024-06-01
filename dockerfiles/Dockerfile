FROM ubuntu:22.04
WORKDIR /
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y
RUN apt-get install build-essential python3-pip python3-dev python3-venv libgtk2.0-dev -y
RUN pip3 install uwsgi flask flask_restful flask-cors pymongo requests boto3 requests celery celery[redis] pika redis
RUN pip3 install flask_jwt_extended uwsgi pyyaml warrant scanpy leidenalg matplotlib
RUN pip3 install opencv-python scipy paramiko pymysql cellpose scikit-image