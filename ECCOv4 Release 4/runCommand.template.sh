export MYQUEUE="--partition standard --qos standard"
export RUNDIR="___RUNDIR___"
myrun=___RUNNAME___
##myargs="--time=12:00:00 --nodes=1 --tasks-per-node=96" #1 node (96 core) mode
##myargs="--time=1-00:00:00 --nodes=2 --tasks-per-node=96" #2 node (192 core) mode
myargs="--time=04:00:00 --nodes=3 --tasks-per-node=121" #3 node (363 core) mode

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ^
# REMEMBER: Update SIZE.h and data.exch2 accordingly |
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ |

sbatch --job-name=$myrun ${myargs//\"} ${MYQUEUE//\"}  ___RUNSUBMIT___

