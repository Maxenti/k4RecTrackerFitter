source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh -r 2025-08-25
mkdir -p /tmp/k4rt_test && cd /tmp/k4rt_test
xrdcp -f root://eosuser.cern.ch//eos/user/c/cglenn/k4rt_dch/k4rt_runtime.tar.gz .
tar xzf k4rt_runtime.tar.gz

which python3
export GAUDI_PLUGIN_PATH=$PWD:$PWD/Tracking:$PWD/DCHdigi
export GAUDI_CONFDB_PATH=$PWD

python3 -c "import TrackingConf; from Configurables import GGTF_tracking; print(GGTF_tracking)"
