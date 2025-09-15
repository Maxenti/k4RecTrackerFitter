ddsim --steeringFile /afs/cern.ch/user/c/cglenn/FCCWork/k4RecTracker/Tracking/test/testTrackFinder/SteeringFile_IDEA_o1_v03DCH.py \
      --compactFile  /eos/user/c/cglenn/FCCWork/GithubRepos/k4geoMax/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03CF_2umAu.xml \
      -G --gun.distribution uniform --gun.particle mu- \
      --random.seed 42 \
      --numberOfEvents 1000 \
      --outputFile /afs/cern.ch/user/c/cglenn/FCCWork/k4RecTracker/Tracking/test/testTrackFinder/out_sim_edm4hep.root 