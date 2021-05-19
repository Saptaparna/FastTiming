#!/bin/tcsh
echo $PWD

echo ${1}
echo ${2}
echo ${3}

cat>Job_${1}_${2}_${3}.csh<<EOF
#!/bin/tcsh
#xrdcp -s root://cmseos.fnal.gov//store/user/sapta/CMSSW_11_1_0_pre7.tar.gz .
tar -xf CMSSW_11_1_0_pre7.tar.gz
source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH slc7_amd64_gcc820 
cd CMSSW_11_1_0_pre7/src
scramv1 b ProjectRename
cmsenv
setenv SCRAM_ARCH "slc7_amd64_gcc820";
scramv1 b
cd -
cmsRun  step1_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg.py >& step1_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}.txt
EOF

chmod 775 Job_${1}_${2}_${3}.csh

cat>condor_${1}_${2}_${3}.jdl<<EOF
universe = vanilla
Executable = Job_${1}_${2}_${3}.csh
Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )
request_disk = 10000000
request_memory = 4000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = Job_${1}_${2}_${3}.csh, step1_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg.py, CMSSW_11_1_0_pre7.tar.gz
notification = Never
Output = CondorJobs/STDOUT_${1}_${2}_${3}_PU.stdout
Error = CondorJobs/STDERR_${1}_${2}_${3}_PU.stderr
Log = CondorJobs/LOG_${1}_${2}_${3}_PU.log
#x509userproxy = ${X509_USER_PROXY}
Queue 1
EOF

condor_submit condor_${1}_${2}_${3}.jdl
