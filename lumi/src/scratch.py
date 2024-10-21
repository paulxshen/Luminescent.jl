import stltovoxel
stltovoxel.convert_files(['in1.stl', 'in.stl'], 'temp/output.png',
                         colors=[(255, 0, 0), (0, 0, 255)], voxel_size=0.01, pad=0)
