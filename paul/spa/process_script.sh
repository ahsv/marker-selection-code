mkdir raw-data
mkdir split-data

for i in {1..10}; 
do  
    cp ${i}/ovrResults-${i}.dat raw-data/; 
    cp ${i}/ovrNumGenes-${i}.dat raw-data/; 
done

cd raw-data

for i in {1..10}; 
do 
    awk -f rowAves.awk ovrNumGenes-${i}.dat > aveGenes-${i}.dat; 
    paste ovrResults-${i}.dat aveGenes-${i}.dat > full-${i}.dat; 
    split -l 11 --numeric-suffixes full-${i}.dat ${i}-; 
    mv ${i}-* ../split-data; 
done

for i in {00..14}; 
do 
    for j in {1..10}; 
    do 
        cat ../split-data/${j}-${i} >> ../lamb-${i}.dat; 
    done; 
    awk '{if ($1 !~ /#/) {print;}}' ../lamb-${i}.dat > lamb-${i}.tmp; 
    mv lamb-${i}.tmp ../lamb-${i}.dat; 
done

cd ..

for i in {00..14};
do
    awk -f procData.awk lamb-${i}.dat >> procData.dat;
done

sed -n '1~3p' procData.dat > effMarks.dat
sed -n '2~3p' procData.dat > maxCorr.dat
sed -n '3~3p' procData.dat > minMarks.dat
