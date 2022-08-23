# a script to check if jobs which are run with submitjobs are completed; hangs on until then or time limit
# takes job id as the only argument

import time
import subprocess
import sys

job_id=sys.argv[1]
t = time.time()

def check_for_completion(job_id):
    jobs=subprocess.check_output(['squeue','-u','groudko']).decode()
    if not job_id in str(jobs):
        return True

while True:
    # Break if this takes more than some_limit
    # 10 hours = 360000 s
    # 48 hours = 1728000 s
    if time.time() - t > 1728000:
        print('time limit (sth could have gone wrong with some samples)',file=sys.stderr)
        subprocess.run(['squeue','-u','groudko','-j'+job_id])
        break
    # Check if the jobs are done. This could be done by
    # grep'ing squeue for your username and some tags
    # that you name your jobs
    if check_for_completion(job_id):
        break
    # Sleep for a while depending on the estimated completion time of the jobs
    time.sleep(60)
