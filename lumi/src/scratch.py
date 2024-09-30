import itertools
from pprint import pprint
import luminescent as lumi
# lumi.make_training_movie(name="mode_converter")
# lumi.make_simulation_movie(name="mode_converter")
name = "mode_converter"
sol = lumi.load_solution(name=name)
sol["optimized_component"].show()
