# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 13:25:21 2021

@author: Owner
"""

#https://stackoverflow.com/questions/52057540/redirecting-pythons-console-output-to-

import osparc
from osparc.models import File, Solver, Job, JobStatus, JobInputs, JobOutputs
from osparc.api import FilesApi, SolversApi
import os
import time
from pathlib import Path
import sys
from page_setup import get_user_ID

os.environ["OSPARC_API_KEY"] = "7d382cb1-0779-562b-8b15-75c9ce7a7178"
os.environ["OSPARC_API_SECRET"] = "cf80455e-deec-5bfc-bbd2-3e87bc8b3ab2"

try: cfg = osparc.Configuration(
              username=os.environ["OSPARC_API_KEY"],
                password=os.environ["OSPARC_API_SECRET"])
except: 
    print("missing os.environ[OSPARC_API_KEY, OSPARC_API_SECRET]")
    cfg = osparc.Configuration(username="*",password='*')


def upload_files(cfg, axons, session, user):
    with osparc.ApiClient(cfg) as api_client:
        files_api = FilesApi(api_client)
        
        #session from input
        #ID  = 1
        path = '../data/u/{}/{}/'.format(user,session)
        
        input_file_1: File = files_api.upload_file(file= path + "array.json")
        input_file_2: File = files_api.upload_file(file=path + "nerve.xml")
        input_file_3: File = files_api.upload_file(file = path + "nerve.json" )
        input_file_4 : File = files_api.upload_file(file = "../data/share/axon/{}.mat".format(axons))
        print('Files uploaded : ', files_api.list_files())

        
        """
        input_file_1 : File = files_api.get_file('09c57c07-da11-3bd0-b554-5705d1ea3055')
        input_file_2 : File = files_api.get_file('5fb397b5-5f21-3ab5-9d09-ad41032793e2')
        input_file_3 : File = files_api.get_file('53949ed5-1305-38c1-b640-d6e8f662d2f8')
        input_file_4 : File = files_api.get_file("e2f44558-ec08-3470-93d3-64c08bc9aa69")
        """
        
        
        return input_file_1, input_file_2, input_file_3, input_file_4
    
#"simcore/services/comp/nerve-mesh", "1.0.2"
def run_node(cfg, inputs, name, version):
    with osparc.ApiClient(cfg) as api_client:

        solvers_api = SolversApi(api_client)

        solver: Solver = solvers_api.get_solver_release(
        name, version)  
        print('Set version and name of solver', solver.id, solver.version)
         
        keys = ['input_' + str(i + 1) for i in range(len(inputs))]        
        if len(keys) == len(inputs):
            pass
        else:
            print('cannot zip')
            sys.exit()        
        job_dict = {}
        for k, v in zip(keys, inputs):
            job_dict.update({k : v})     
        
     
        job: Job = solvers_api.create_job(
        solver.id,
        solver.version,
        JobInputs(
               job_dict
             )
           )
        print('Created Job', job.id)
        
        print('Starting job and printing status')
        status: JobStatus = solvers_api.start_job(solver.id, solver.version, job.id)    
        while not status.stopped_at:
            time.sleep(5)
            status = solvers_api.inspect_job(solver.id, solver.version, job.id)
            print("Solver progress", f"{status.progress}/100", flush=True)
                    
        #Outputs
        outputs: JobOutputs = solvers_api.get_job_outputs(solver.id, solver.version, job.id)
                  
        #Print outputs
        print(f"Job {outputs.job_id} got these results:")
        for output_name, result in outputs.results.items():
            print(output_name, "=", result)
        
        files_api = FilesApi(api_client)
        results_file: File = outputs.results["output_1"]
        download_path: str = files_api.download_file(file_id=results_file.id)
        #print(Path(download_path).read_text())
            
        return results_file, download_path 
 
    
#testing the api; comment later
"""
#step 1
#upload files
input_file_1, input_file_2, input_file_3, input_file_4 = upload_files(cfg)

#step 2 
#run nerve_mesh

nm_results, nm_download_path = run_node(cfg, \
         [input_file_1, input_file_2, input_file_3],\
         "simcore/services/comp/nerve-mesh", \
         "1.0.2")

#run axon-population
ap_results, ap_download_path = run_node(cfg, \
         [input_file_4, input_file_2, input_file_3],\
         "simcore/services/comp/axon-population", \
         "1.0.2")

#run
es_results, es_download_path = run_node(cfg, \
         [nm_results],\
         "simcore/services/comp/eidors-solver", \
         "2.0.2")

#nerve recording
nr_results, nr_download_path = run_node(cfg, \
         [ap_results, es_results],\
         "simcore/services/comp/nerve-recording", \
         "1.0.1")

#es results and nr_results must be vis

"""

    
    
    
  

