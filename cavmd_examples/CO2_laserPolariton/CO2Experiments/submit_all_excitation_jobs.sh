#/bin/bash

for path in $(ls -d exciteCO2*)
do
    echo "Entering $path"
    cd $path

    ./submit_jobs.sh
    
    cd ..
done


wait
