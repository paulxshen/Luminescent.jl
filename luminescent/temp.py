# import skrf as rf
# a = rf.Network('S.s3p')
# s = a.s
# f = a.f
# s = s[:, :2, :2]
# f = rf.Frequency.from_f(f, unit='GHz')
# a = rf.Network(frequency=f, s=s)
# a.write_touchstone('S.s2p')

import luminescent as lumi
# lumi.make_simulation_movie(     "build/precompile_execution/tiny_2_float32_CUDA")
# lumi.finetune("build/precompile_execution/back_float32", 1, framerate=20)
lumi.make_simulation_movie("build/precompile_execution/back_float32")
