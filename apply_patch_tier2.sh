#!/bin/bash

sed -i 's/images.dodas.infn.it/unpacked.cern.ch/g' /usr/local/lib/python3.9/site-packages/dask_remote_jobqueue/templates/scheduler.sub
sed -i 's/{{ selected_sitename }}/requirements = ( SiteName == "T2_LNL_PD_Af20" )/g' /usr/local/lib/python3.9/site-packages/dask_remote_jobqueue/templates/scheduler.sub
sed -i 's/ON_EXIT/ON_EXIT_OR_EVICT/g' /usr/local/lib/python3.9/site-packages/dask_remote_jobqueue/templates/start_scheduler.py
sed -i 's/condor_rm/condor_suspend/g' /usr/local/lib/python3.9/site-packages/dask_remote_jobqueue/templates/job_rm.sh
