#!/bin/bash

chmod 600 .globus/usercert.pem
chmod 600 .globus/userkey.pem

voms-proxy-init --voms cms --key .globus/userkey.pem --cert .globus/usercert.pem
