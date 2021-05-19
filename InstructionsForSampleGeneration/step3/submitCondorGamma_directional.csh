#!/bin/tcsh
echo $PWD

echo ${1}
echo ${2}
echo ${3}

cat>Job_${1}_${2}_${3}_directional.csh<<EOF
#!/bin/tcsh
#xrdcp -s root://cmseos.fnal.gov//store/user/sapta/CondorTar/CMSSW_11_1_0_pre7.tar.gz .
tar -xf CMSSW_11_1_0_pre7.tar.gz
source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH slc7_amd64_gcc820 
cd CMSSW_11_1_0_pre7/src
scramv1 b ProjectRename
cmsenv
setenv SCRAM_ARCH "slc7_amd64_gcc820";
scramv1 b
cd -
cmsRun  step3_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_directional.py >& step3_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_directional.txt
EOF

chmod 775 Job_${1}_${2}_${3}_directional.csh

cat>condor_${1}_${2}_${3}.jdl<<EOF
universe = vanilla
Executable = Job_${1}_${2}_${3}_directional.csh
Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )
request_disk = 10000000
request_memory = 8000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = Job_${1}_${2}_${3}_directional.csh, step3_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_directional.py, CMSSW_11_1_0_pre7.tar.gz
notification = Never
Output = CondorJobs/STDOUT_${1}_${2}_${3}_directional_PU.stdout
Error = CondorJobs/STDERR_${1}_${2}_${3}_directional_PU.stderr
Log = CondorJobs/LOG_${1}_${2}_${3}_directional_PU.log
x509userproxy = ${X509_USER_PROXY}
Queue 1
EOF

condor_submit condor_${1}_${2}_${3}.jdl
