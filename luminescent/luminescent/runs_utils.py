# import dill
from pprint import pprint
import os
import json
import numpy as np
from .pic.inverse_design import *
from .pic.sparams import *
from .utils import *
from .layers import *
from .constants import *


def lastrun(name="", wd=os.path.join(os.getcwd(), "runs"), study="",  **kwargs):
    if name:
        return os.path.join(wd, name)
    l = [os.path.join(wd, x) for x in os.listdir(wd)]
    l = [x for x in l if os.path.isfile(os.path.join(x, "problem.json"))]
    l = sorted(l, key=lambda x: os.path.getmtime(x), reverse=True)
    if study:
        for x in l:
            try:
                if json.loads(open(os.path.join(x, "problem.json"), "rb").read())["study"] == study:
                    return x
            except:
                pass
    return l[0]


def load_prob(path):
    path = os.path.abspath(path)
    print(f"loading problem from {path}")
    return json.loads(open(os.path.join(path, "problem.json"), "rb").read())


def make_simulation_movie(framerate=30, **kwargs):
    1
    # path = lastrun(**kwargs)

    # f = os.path.join(path, "geometry")
    # imgs = sorted(os.listdir(f), key=lambda x: float(x[0:-4]))
    # frame = cv2.imread(os.path.join(f, imgs[0]))
    # height, width, layers = frame.shape

    # video = cv2.VideoWriter(os.path.join(
    #     path, "simulation_video.mp4"), 0x7634706d, framerate, (width, height))

    # for img in imgs:
    #     video.write(cv2.imread(os.path.join(f, img)))

    # cv2.destroyAllWindows()
    # video.release()


def make_training_movie(framerate=2, **kwargs):
    1
    # path = lastrun(**kwargs)

    # f = os.path.join(path, "checkpoints")
    # ckpts = sorted(os.listdir(f))
    # frame = cv2.imread(os.path.join(f, ckpts[0], "run_1.png"))
    # height, width, layers = frame.shape

    # video = cv2.VideoWriter(os.path.join(
    #     path, "training_video.mp4"), 0x7634706d, framerate, (width, height))

    # for ckpt in ckpts:
    #     video.write(cv2.imread(os.path.join(f, ckpt, "run_1.png")))

    # cv2.destroyAllWindows()
    # video.release()
