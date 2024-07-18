import os
import subprocess
import time
import json
import bson
import numpy as np
from .inverse_design import *
from .sparams import *
from .utils import *
from .layers import *
from .constants import *


from subprocess import Popen, PIPE


def solve(prob, ):
    c0 = prob["component"]
    del prob["component"]
    bson_data = bson.dumps(prob)
    prob["component"] = c0

    path = prob["path"]
    if not os.path.exists(path):
        os.makedirs(path)
    print(f"""
          using simulation folder {path}
          started julia process
          """)
    prob_path = os.path.join(path, "prob.bson")
    with open(prob_path, "wb") as f:
        # Write the BSON data to the file
        f.write(bson_data)

    start_time = time.time()

    cmd = [f"julia", os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "run.jl"), prob_path,]
    with Popen(cmd,  stdout=PIPE, ) as p:
        for line in p.stdout:
            print(line, flush=True)
    # exit_code = p.poll()
    # subprocess.run()
    # print(f"julia simulation took {time.time()-start_time} seconds")
    print(f"images and results saved in {path}")
    sol = load_solution(path)
    if prob["study"] == "inverse_design":
        c = apply_design(c0,  sol)
        sol["before"]["component"] = c0
        sol["after"]["component"] = c
    return sol


def load_sparams(sparams):
    if "re" in list(sparams.values())[0]:
        return {k: v["re"]+1j*v["im"] for k, v in sparams.items()}
    return {wl: {k: (v["re"]+1j*v["im"])
                 for k, v in d.items()} for wl, d in sparams.items()}


def load_solution(path=None):
    if path is None:
        path = sorted(os.listdir(PATH))[-1]
        path = os.path.join(PATH, path)
    print(f"loading solution from {path}")
    prob = bson.loads(open(os.path.join(path, "prob.bson"), "rb").read())
    p = os.path.join(path, "sol.json")
    # sol = bson.loads(p, "rb").read())["sol"]
    sol = json.loads(open(p).read())
    if prob["study"] == "sparams":
        sol["sparams"] = load_sparams(sol["sparams"])
    elif prob["study"] == "inverse_design":
        for k in ["before", "after"]:
            sol[k]["sparams"] = load_sparams(sol[k]["sparams"])
        for k in sol["after"]:
            sol[k] = sol["after"][k]
    # pic2gds(os.path.join(path, f"design{i+1}.png"), sol["dx"])

        sol["designs"] = [np.array(d) for d in sol["designs"]]
        pass
    return sol


def write_sparams(*args, **kwargs):
    return solve(sparams_problem(*args, **kwargs))
