#!/bin/tcsh
echo $PWD

echo ${1}
echo ${2}
echo ${3}

cat>Job_${1}_${2}_${3}_directional.csh<<EOF
#!/bin/tcsh
tar -xzf CMSSW_9_1_0_pre3.tar.gz
cd CMSSW_9_1_0_pre3
scram b ProjectRename
source /cvmfs/cms.cern.ch/cmsset_default.csh
cmsenv
cd -
cmsRun  step3_SinglekLongGun_BatchSubmit_NEVENT${1}_NPART${2}_PT${3}_cfg_directional_PU.py >& step3_SinglekLongGun_BatchSubmit_NEVENT${1}_NPART${2}_PT${3}_directional_PU.txt
EOF

chmod 775 $PWD/Job_${1}_${2}_${3}_directional.csh

cat>condor_${1}_${2}_${3}.jdl<<EOF
universe = vanilla
Executable = Job_${1}_${2}_${3}_directional.csh
Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )
request_disk = 10000000
request_memory = 8100
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = /uscms_data/d2/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/Restart/CMSSW_9_1_0_pre3/src/CMSSW_9_1_0_pre3.tar.gz, /uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/Restart/CMSSW_9_1_0_pre3/src/kLongTime_step3_PU/Job_${1}_${2}_${3}_directional.csh, /uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/Restart/CMSSW_9_1_0_pre3/src/kLongTime_step3_PU/step3_SinglekLongGun_BatchSubmit_NEVENT${1}_NPART${2}_PT${3}_cfg_directional_PU.py
notification = Never
Output = CondorJobs/STDOUT_${1}_${2}_${3}_directional_PU.stdout
Error = CondorJobs/STDERR_${1}_${2}_${3}_directional_PU.stderr
Log = CondorJobs/LOG_${1}_${2}_${3}_directional_PU.log
x509userproxy = ${X509_USER_PROXY}
Queue 1
EOF

condor_submit condor_${1}_${2}_${3}.jdl
