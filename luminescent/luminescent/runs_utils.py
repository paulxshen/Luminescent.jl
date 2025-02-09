# import dill
from pprint import pprint
import os
import json
import cv2
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


def load_problem(path):
    path = os.path.abspath(path)
    print(f"loading problem from {path}")
    return json.loads(open(os.path.join(path, "problem.json"), "rb").read())


def make_simulation_movie(path, framerate=None, **kwargs):
    TEMP = os.path.join(path, "temp")
    FIELDS = os.path.join(TEMP, "fields")
    FRAMES = os.path.join(TEMP, "frames")
    MOVIE = os.path.join(path, "simulation_video.mp4")
    shutil.rmtree(FRAMES, ignore_errors=True)
    os.makedirs(FRAMES, exist_ok=True)
    shutil.rmtree(MOVIE, ignore_errors=True)

    if framerate is None:
        framerate = load_problem(path)["framerate"]

    g = np.load(os.path.join(path, "temp", 'g.npy')).T
    gmax = np.max(np.abs(g))

    dir = os.path.join(path, "temp", 'fields')
    umax = 0
    fns = sorted(os.listdir(dir), key=lambda x: float(x[0:-4]))
    for fn in fns:
        a = np.load(os.path.join(dir, fn))
        v = np.max(np.abs(a))
        if umax < v:
            umax = v

    for fn in fns:
        name = fn[0:-4]
        a = np.load(os.path.join(FIELDS, fn)).T
        fig, axs = plt.subplots(1, 2)
        axs[1].imshow(a, cmap='seismic', origin='lower',
                      vmin=-umax, vmax=umax)
        axs[0].imshow(-g, cmap='gray',
                      origin='lower', vmin=-gmax, vmax=0)
        plt.savefig(os.path.join(FRAMES, f"{name}.png"))
        plt.close(fig)

    fns = sorted(os.listdir(FRAMES), key=lambda x: float(x[0:-4]))
    frame = height = width = layers = video = None
    for i, fn in enumerate(fns):
        frame = cv2.imread(os.path.join(FRAMES, fn))
        if i == 0:
            height, width, layers = frame.shape
            video = cv2.VideoWriter(
                MOVIE, 0x7634706d, framerate, (width, height))

        video.write(frame)

    cv2.destroyAllWindows()
    video.release()


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
