import dill
from PIL import Image
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


def solve(prob, dev=False, run=True):
    if "dev" in prob:
        dev = prob["dev"]
    # c0 = prob["component"]
    # del prob["component"]
    bson_data = bson.dumps(prob)
    # prob["component"] = c0

    path = prob["path"]
    if not os.path.exists(path):
        os.makedirs(path)
    print(f"""
          using simulation folder {path}
          started julia process
          compiling julia code...
          """)
    prob_path = os.path.join(path, "prob.bson")
    with open(prob_path, "wb") as f:
        # Write the BSON data to the file
        f.write(bson_data)

    if not run:
        return
    start_time = time.time()

    def run(cmd):

        proc = Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # proc.wait()
        with proc:
            for line in proc.stdout:

                print(str(line.decode().strip()), flush=True)
            err_message = proc.stderr.read().decode()
            print(err_message)
    cmd_dev = [f"julia", os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "run.jl"), path,]
    cmd = ["lumi", path]
    if dev:
        print(" ".join(cmd_dev))
        run(cmd_dev)
    else:
        try:
            run(cmd)
        except Exception as e:
            print(e)
            run(cmd_dev)

    # with Popen(cmd,  stdout=PIPE, stderr=PIPE) as p:
    #     if p.stderr is not None:
    #         for line in p.stderr:
    #             print(line, flush=True)
    # exit_code = p.poll()
    # subprocess.run()
    # print(f"julia simulation took {time.time()-start_time} seconds")
    print(f"images and results saved in {path}")
    sol = load_solution(path=path)
    if prob["study"] == "inverse_design":
        c = apply_design(c0,  sol)
        sol["optimized_component"] = c
    #     sol["before"]["component"] = c0
    #     sol["after"]["component"] = c
    return sol


def get_path(path=None, study=""):
    if path is None:
        l = sorted(os.listdir(RUNS_PATH), reverse=True)
        if study:
            for p in l:
                try:
                    s = json.loads(open(os.path.join(RUNS_PATH, p, "sol.json")).read())[
                        "study"]
                    if s == study:
                        path = p
                        break
                except:
                    pass
        else:
            path = l[0]
        path = os.path.join(RUNS_PATH, path)
        return path
    if os.path.isdir(path):
        return path

    if os.path.isdir(os.path.join(RUNS_PATH, path)):
        return path


def load_component(path=None, study=None):
    return dill.load(os.path.join(get_path(path), "comp.pk"))


def finetune(iters, path=None,):
    if path is None:
        path = get_path(study="inverse_design")

    prob = bson.loads(open(os.path.join(path, "prob.bson"), "rb").read())
    prob["iters"] = iters
    # with open(path, "wb") as f:
    #     f.write(bson.dumps(prob))
    solve(prob)


def load_sparams(sparams):
    if "re" in list(sparams.values())[0]:
        return {k: v["re"]+1j*v["im"] for k, v in sparams.items()}
    return {wl: {k: (v["re"]+1j*v["im"])
                 for k, v in d.items()} for wl, d in sparams.items()}


def load_solution(path=None, study="",):
    if path is None:
        path = get_path(study=study)
    print(f"loading solution from {path}")
    prob = bson.loads(open(os.path.join(path, "prob.bson"), "rb").read())
    p = os.path.join(path, "sol.json")
    # sol = bson.loads(p, "rb").read())["sol"]
    sol = json.loads(open(p).read())
    sol["sparams"] = load_sparams(sol["sparams"])
    if prob["study"] == "sparams":
        pass
    elif prob["study"] == "inverse_design":
        # for k in ["before", "after"]:
        #     sol[k]["sparams"] = load_sparams(sol[k]["sparams"])
        # for k in sol["after"]:
        #     sol[k] = sol["after"][k]

        l = [np.array(d) for d in sol["optimized_designs"]]
        sol["optimized_designs"] = l
        for i, a in enumerate(l):
            Image.fromarray(np.uint8(a * 255),
                            'L').save(os.path.join(path, f"design{i+1}.png"))
            pic2gds(os.path.join(path, f"design{i+1}.png"), sol["dx"])
        pass
    return sol


def write_sparams(*args, run=True, **kwargs):
    return solve(sparams_problem(*args, **kwargs), run=run)
