#!/bin/bash
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Li 1 2 6 6 1.5 1' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Li 0 2 6 6 1.5 1' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Li 1 2 6 6 6.0 1' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Li 0 2 6 6 6.0 1' &  
wait
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Otsu 1 2 6 6 1.5 1' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Otsu 0 2 6 6 1.5 1' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Otsu 1 2 6 6 6.0 1' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Otsu 0 2 6 6 6.0 1' &  
wait

nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Li 1 2 6 6 1.5 0' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Li 0 2 6 6 1.5 0' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Li 1 2 6 6 6.0 0' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Li 0 2 6 6 6.0 0' &  
wait
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Otsu 1 2 6 6 1.5 0' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Otsu 0 2 6 6 1.5 0' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Otsu 1 2 6 6 6.0 0' &  
nohup nice -n 5 ImageJ-linux64 --allow-multiple --headless -macro ./process_embryo.ijm  'Otsu 0 2 6 6 6.0 0' &  
exit
