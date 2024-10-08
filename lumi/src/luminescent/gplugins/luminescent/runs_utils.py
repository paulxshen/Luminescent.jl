# import dill
import cv2
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
from .sol import *


def lastrun(wd=os.path.join(os.getcwd(), "runs"), name="", study="",  **kwargs):
    if name:
        return os.path.join(wd, name)
    l = [os.path.join(wd, x) for x in os.listdir(wd)]
    l = [x for x in l if os.path.isfile(os.path.join(x, "problem.bson"))]
    l = sorted(l, key=lambda x: os.path.getmtime(x), reverse=True)
    if study:
        for x in l:
            try:
                if bson.loads(open(os.path.join(x, "problem.bson"), "rb").read())["study"] == study:
                    return x
            except:
                pass
    return l[0]


def finetune(iters, **kwargs):
    path = lastrun(study="inverse_design", **kwargs)

    prob = bson.loads(open(os.path.join(path, "problem.bson"), "rb").read())
    prob["iters"] = iters
    # prob["eta"] = eta
    prob["restart"] = False
    prob = {**prob, **kwargs}
    return solve(prob)


def load_problem(**kwargs):
    path = lastrun(**kwargs)
    print(f"loading problem from {path}")
    return bson.loads(open(os.path.join(path, "problem.bson"), "rb").read())


def load_solution(**kwargs):
    path = lastrun(**kwargs)
    print(f"loading solution from {path}")
    prob = bson.loads(open(os.path.join(path, "problem.bson"), "rb").read())
    p = os.path.join(path, "solution.json")
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
            Image.fromarray(np.flip(np.uint8((1-a)) * 255, 0),
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


def make_simulation_movie(framerate=30, **kwargs):
    path = lastrun(**kwargs)

    f = os.path.join(path, "temp")
    imgs = sorted(os.listdir(f), key=lambda x: float(x[0:-4]))
    frame = cv2.imread(os.path.join(f, imgs[0]))
    height, width, layers = frame.shape

    video = cv2.VideoWriter(os.path.join(
        path, "simulation_video.mp4"), 0x7634706d, framerate, (width, height))

    for img in imgs:
        video.write(cv2.imread(os.path.join(f, img)))

    cv2.destroyAllWindows()
    video.release()


def make_training_movie(framerate=2, **kwargs):
    path = lastrun(**kwargs)

    f = os.path.join(path, "checkpoints")
    ckpts = sorted(os.listdir(f))
    frame = cv2.imread(os.path.join(f, ckpts[0], "run_1.png"))
    height, width, layers = frame.shape

    video = cv2.VideoWriter(os.path.join(
        path, "training_video.mp4"), 0x7634706d, framerate, (width, height))

    for ckpt in ckpts:
        video.write(cv2.imread(os.path.join(f, ckpt, "run_1.png")))

    cv2.destroyAllWindows()
    video.release()


def write_sparams(*args, run=True, **kwargs):
    return solve(sparams_problem(*args, **kwargs), run=run)
