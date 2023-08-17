chmod +x job_submit.sh
chmod +x job_rm.sh
chmod +x entrypoint.sh


while true; do
    curl -d grant_type=urn:ietf:params:oauth:grant-type:token-exchange \
        -u $IAM_CLIENT_ID:$IAM_CLIENT_SECRET \
        -d audience="https://wlcg.cern.ch/jwt/v1/any" \
        -d subject_token=`cat token` \
        -d scope="openid profile wlcg wlcg.groups" \
        ${IAM_SERVER}/token \
        | tee /tmp/response | jq .access_token |  tr -d '"' |  tr -d '\n'> /tmp/token_tmp \
    && cp /tmp/token_tmp token
    sleep 72000
done &


python3 start_scheduler.py
