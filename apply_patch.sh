#!/bin/bash

cp ./patches/scheduler.sub /usr/local/lib/python3.9/site-packages/dask_remote_jobqueue/templates/scheduler.sub 
cp ./patches/scheduler.sh /usr/local/lib/python3.9/site-packages/dask_remote_jobqueue/templates/scheduler.sh
