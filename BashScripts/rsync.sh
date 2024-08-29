#!/bin/bash

ssh curie "cd /home/data/PD/Simulations/NewRuns && find" | grep -e 'results/\S*.jld2' -e 'JobData/\S*.jld2' | rsync -ravz --files-from=- curie:/home/data/PD/Simulations/NewRuns ~/Projects/TSF_Notebooks//NewRuns
