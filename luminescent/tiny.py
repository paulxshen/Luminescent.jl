import os
import luminescent as lumi
import gdsfactory as gf

# dir = "runs"
dir = os.path.join("build", "precompile_execution")

c = gf.components.straight(length=.1, width=0.5,)
wavelengths = 1.55

path = os.path.join(dir, "tiny")
lumi.make_pic_sim_prob(path, c, wavelengths=wavelengths, keys=[
                       "2,1"], nres=15, approx_2D_mode="TE")
path = os.path.join(dir, "tinycu")
lumi.make_pic_sim_prob(path, c, wavelengths=wavelengths, keys=[
                       "2,1"], nres=15, approx_2D_mode="TE", gpu="CUDA")

path = os.path.join(dir, "tiny3")
lumi.make_pic_sim_prob(path, c, wavelengths=wavelengths, keys=[
                       "2,1"],                       nres=15, )
path = os.path.join(dir, "tiny3cu")
lumi.make_pic_sim_prob(path, c, wavelengths=wavelengths, keys=[
                       "2,1"],                       nres=15, gpu="CUDA")


path = os.path.join(dir, "back")
# lumi.load_res()
# raise NotImplementedError("This is a stub")
c = lumi.mimo(west=1, east=1, l=.5, w=.5,  wwg=.5)
targets = {"tparams": {
    1.55: {
        "2,1": 1.0
    }}}
lumi.make_pic_inv_prob(
    path,  c, targets,
    lvoid=0.2, iters=2, nres=15,
    approx_2D_mode="TE")


# lumi.finetune(os.path.join("runs", "back"), iters=2)
