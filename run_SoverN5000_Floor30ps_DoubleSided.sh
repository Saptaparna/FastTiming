echo "Now running SoverN5000 Floor30ps Double Sided"

mkdir OutputFiles/SoverN5000_Floor30ps_DoubleSided 

for i in $(seq 0.1 0.1 0.9)
do
  for j in $(seq 0.1 0.1 0.9)
  do
  echo ${i}
  echo ${j}
  nohup root -l -b -q 'TruncatedMean_DoubleSided.cc++("InputFiles/HGCTiming_kLong_Pt10_SoverN5000ps_Floor30ps", "OutputFiles/SoverN5000_Floor30ps_DoubleSided/output_HGCTiming_kLong_Pt10_SoverN5000ps_Floor30ps", 5.0, '${i}', '${j}')' > OutputFiles/SoverN5000_Floor30ps_DoubleSided/output_HGCTiming_kLong_Pt10_SoverN5000ps_Floor30ps_FractionLow_${i}_FractionHigh_${j}.txt
  done
done
