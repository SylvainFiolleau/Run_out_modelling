# This script can run within QGIS python console, but also with your prefered python method if AvaFrame is not installed via QGIS.
# The script runs a bunch of AvaFrame simulations in batch (all those contained in the out_folder and sub scenario folders)
# Remember to update the path in the local_avaframeCfg.ini file of each simulation scenario contained in the out_folder/sub_folders if files/folders are moved around!
# Author: Francois Noel, 2024. Open source CC BY 4.0 (free to use, modify, sell if citing original author)

import os


def run_avaframe(out_folder, ava_folder):
    print('Start simulation')

    # if True: #simple switch to only print the name of the remaining sites (False) or to also simulate them (True)

    path_site = out_folder

    # read the custom AvaFrame settings file and write it to the Avaframe com1DFA folder
    with open(os.path.join(path_site, 'local_avaframeCfg.ini'), 'r') as fp:
        data = fp.readlines()  # read a list of lines into data
    with open(os.path.join(ava_folder, 'local_avaframeCfg.ini'), 'w') as fp:
        fp.writelines(data)  # write the AvaFrame parameter file into the com1DFA folder

    # read the custom AvaFrame parameter file and write it to the Avaframe com1DFA folder
    with open(os.path.join(path_site, 'local_com1DFACfg.ini'), 'r') as fp:
        data = fp.readlines()  # read a list of lines into data
    with open(os.path.join(ava_folder, 'com1DFA', 'local_com1DFACfg.ini'), 'w') as fp:
        fp.writelines(data)  # write the AvaFrame parameter file into the com1DFA folder

    # run Avaframe
    dfa = "python " + os.path.join(ava_folder, 'runCom1DFA.py')
    os.system(dfa)

    print('Done running the AvaFrame simulation!\n')
