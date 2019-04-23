import read_MLFMA_mesh_part1
import mesh_functions_seb
import os

simuDirName = "."
inputDirName = "./inputParams"

files = os.listdir("./mesh/")
# for matFileName in files:
#     os.remove("./mesh/" + matFileName)
read_MLFMA_mesh_part1.run(simuDirName, inputDirName)
mesh_functions_seb.run(simuDirName)
