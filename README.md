# AtlasWorkers

AtlasWorkers is celery-based distributed computation module for AtlasWebService. However, it is a general computation module for other tools with standard interface by using RabbitMQ, Redis.

### Recommended Usage

##### Install and execution

Before running, `config.yml` file is required to be ready with AWS credentials as well as other options. `config_template.yml` can be copied and furnished with required parameters.

Using docker-compose (or docker compose in recent docker versions). Exmemplary docker-compose.yml is in `dockerfiles/docker-compose.yml`
```
    $ docker-compose run <worker-service-name>
```

For service
```
    $ docker-compose up -d <worker-service-name>
```

##### API endpoint

- [POST] /api/v1/task - Send a worker request with following arguments in JSON format
    - task : task name (e.g. gene.compute_qc)
    - queue : queue name (e.g. atxcloud_gene)
    - args : list of arguments for the task
    - kwargs : keyword arguments for the task

- [GET] /api/v1/task/\<task_id\> - To get an progress or the result of the request. 

### List of current workers

- gene (queue : atxcloud_gene)
    - **compute_qc** , this is used for the AtlasGX gene viewer
        - arguments : [filename, gene_list], e.g. ["data/D172/out/Gene/raw/spatial/genes.h5ad",["Vegfa"]]
        - example payload : 
            ```
                {
                    "queue":"atxcloud_gene",
                    "task":"gene.compute_qc",
                    "args":["data/D172/out/Gene/raw/spatial/genes.h5ad",["Vegfa"]],
                    "kwargs":{}
                }
            ```
        - example response : 
            ```
                {
                    "_id": "218e2d6d-3344-4295-b318-e97c460ea59a",
                    "name": "gene.compute_qc",
                    "args": [
                        "data/D172/out/Gene/raw/spatial/genes.h5ad",
                        [
                            "Vegfa"
                        ]
                    ],
                    "kwargs": {},
                    "queue": "atxcloud_gene",
                    "requested_by": "admin",
                    "requested_at": "2021-11-15T15:24:10.711350"
                }       
            ```
- core (queue : atxcloud_core)
    - This is test purpose worker.

### Changelog

##### 2021-11-15 

- **gene_worker** added with compute_qc task

##### 2021-11-13

- API endpoints are ready in AtlasWebService
- RabbitMQ, Redis deployed in the server (address: api.atlasxomics.com), which are password protected
- Celery added to the docker environment
