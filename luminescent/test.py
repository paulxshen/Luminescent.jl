import os
import luminescent as lumi
import gdsfactory as gf

# dir = "runs"
dir = os.path.join("build", "precompile_execution")
for name in ["tiny", "tinycu", ]:
    path = os.path.join(dir, name)
    lumi.solve(path)
