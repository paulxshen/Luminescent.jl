# import dill
from pprint import pprint
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
    cmd = ["lumi", path]
    run(cmd)

    # with Popen(cmd,  stdout=PIPE, stderr=PIPE) as p:
    #     if p.stderr is not None:
    #         for line in p.stderr:
    #             print(line, flush=True)
    # exit_code = p.poll()
    # subprocess.run()
    # print(f"julia simulation took {time.time()-start_time} seconds")
    print(f"images and results saved in {path}")
    sol = load_solution(path=path)
    return sol


def lastrun(wd=os.path.join(os.getcwd(), "runs"), name="", study="",  **kwargs):
    if name:
        return os.path.join(wd, name)
    l = [os.path.join(wd, x) for x in os.listdir(wd)]
    l = [x for x in l if os.path.isfile(os.path.join(x, "prob.bson"))]
    l = sorted(l, key=lambda x: os.path.getmtime(x), reverse=True)
    if study:
        for x in l:
            try:
                if bson.loads(open(os.path.join(x, "prob.bson"), "rb").read())["study"] == study:
                    return x
            except:
                pass
    return l[0]


def finetune(iters, eta=0.2, **kwargs):
    path = lastrun(study="inverse_design", **kwargs)

    prob = bson.loads(open(os.path.join(path, "prob.bson"), "rb").read())
    prob["iters"] = iters
    prob["eta"] = eta
    prob["restart"] = False
    prob = {**prob, **kwargs}
    return solve(prob)


def load_sparams(sparams):
    if "re" in list(sparams.values())[0]:
        return {k: v["re"]+1j*v["im"] for k, v in sparams.items()}
    return {wl: {k: (v["re"]+1j*v["im"])
                 for k, v in d.items()} for wl, d in sparams.items()}


def load_solution(**kwargs):
    path = lastrun(**kwargs)
    print(f"loading solution from {path}")
    prob = bson.loads(open(os.path.join(path, "prob.bson"), "rb").read())
    p = os.path.join(path, "sol.json")
    # sol = bson.loads(p, "rb").read())["sol"]
    sol = json.loads(open(p).read())
    sol["sparams"] = load_sparams(sol["sparams"])
    sol["component"] = gf.import_gds(os.path.join(path, "component.gds"))
    if prob["study"] == "sparams":
        pass
    elif prob["study"] == "inverse_design":
        l = [np.array(d) for d in sol["optimized_designs"]]
        sol["optimized_designs"] = l
        for i, a in enumerate(l):
            name = f"optimized_design_region_{i+1}.png"
            Image.fromarray(np.uint8((1-a) * 255),
                            'L').save(os.path.join(path, name))
            pic2gds(os.path.join(
                path, name), sol["dx"])
        c = apply_design(sol["component"],  sol)
        # sol["optimized_component"] = copy.deepcopy(c)
        sol["optimized_component"] = c
        c.write_gds(os.path.join(path, "optimized_component.gds"))
    return sol


def show_solution(**kwargs):
    path = lastrun(**kwargs)
    print(f"showing solution from {path}")
    sol = load_solution(**kwargs)
    sol = {k: sol[k] for k in ["path", "sparams", "tparams", ]}
    pprint(sol)

    i = 1
    while True:
        p = os.path.join(path, f"run_{i}.png")
        if os.path.exists(p):
            img = Image.open(p)
            img.show()
            try:
                display(img)
            except:
                pass
            i += 1
        else:
            break


def write_sparams(*args, run=True, **kwargs):
    return solve(sparams_problem(*args, **kwargs), run=run)
