# import dill
from pprint import pprint
import os
import subprocess
from .runs_utils import *
import json
import numpy as np
import requests
from .pic.inverse_design import *
from .pic.sparams import *
from .utils import *
from .layers import *
from .constants import *


from subprocess import Popen, PIPE
URL = "http://127.0.0.1:8975"


def start_fdtd_server(url=URL):
    try:
        r = requests.get(url)
    except:
        print("starting julia fdtd server on localhost...")
        cmd = ["julia", "-e", "using Luminescent; start_fdtd_server()"]
        proc = Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with proc:
            for line in proc.stdout:
                print(str(line.decode().strip()), flush=True)
            err_message = proc.stderr.read().decode()
            print(err_message)


def solve(path):
    prob = load_prob(path)

    print("no fdtd binaries found - starting julia session to compile fdtd code...")
    # prob["action"] = "solve"
    # r = requests.post(f"{url}/local", json=prob)
    1
    # print(f"""
    #       using simulation folder {path}
    #       starting julia process...
    #       """)
    # start_time = time.time()

    def run(cmd):

        proc = Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # proc.wait()
        with proc:
            for line in proc.stdout:

                print(str(line.decode().strip()), flush=True)
            err_message = proc.stderr.read().decode()
            print(err_message)
    env = r'using Pkg;Pkg.activate(raw"c:\Users\pxshe\OneDrive\Desktop\beans\Luminescent.jl\luminescent");'
    env = '0;'
    # cmd = ["lumi", path]
    gpu_backend = prob["gpu_backend"]
    if not gpu_backend:
        cmd = ["julia", "-e", f'{env}using Luminescent;picrun(raw"{path}")']
    else:
        print(f"using {gpu_backend} backend.")
        if gpu_backend == "CUDA":
            cmd = ["julia", "-e",
                   f"{env}using Luminescent,CUDA;@assert CUDA.functional();picrun(\"{path}\";gpuarray=cu)"]
    run(cmd)

    # with Popen(cmd,  stdout=PIPE, stderr=PIPE) as p:
    #     if p.stderr is not None:
    #         for line in p.stderr:
    #             print(line, flush=True)
    # exit_code = p.poll()
    # subprocess.run()
    # print(f"julia simulation took {time.time()-start_time} seconds")
    # print(f"images and results saved in {path}")
    # sol = load_res(path=path)
    # return sol


def load_sparams(sparams):
    if "re" in list(sparams.values())[0]:
        return {k: v["re"]+1j*v["im"] for k, v in sparams.items()}
    return {wl: {k: (v["re"]+1j*v["im"])
                 for k, v in d.items()} for wl, d in sparams.items()}
