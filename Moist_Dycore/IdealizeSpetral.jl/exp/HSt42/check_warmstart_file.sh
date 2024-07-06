#!/bin/bash
if [ -f "RH80_PR5_10day_startfrom_0day_final.dat" ] && [ -f "RH80_PR5_10day_startfrom_2day_final.dat" ] && [ -f "RH80_PR5_10day_startfrom_4day_final.dat" ] && [ -f "RH80_PR5_10day_startfrom_6day_final.dat" ] && [ -f "RH80_PR5_10day_startfrom_8day_final.dat" ]; then
	echo "All files have completed!!!"
	All_file_exist = $true
else
	All_file_exist = $false
fi


while [ "$All_file_exist" ]

julia Run_HS.jl # first run
do 
	if [ -f "warmstart.dat" ]; then
    		echo 'warmstart file exists.'
		### Changing file name
		if [ -f "RH80_PR5_10day_startfrom_0day_final.dat" ] && [ ! -f "RH80_PR5_10day_startfrom_2day_final.dat" ] ; then
			mv "warmstart.dat" "RH80_PR5_10day_startfrom_2day_final.dat"
		
		elif [ -f "RH80_PR5_10day_startfrom_2day_final.dat" ] && [ ! -f "RH80_PR5_10day_startfrom_4day_final.dat" ]; then
                        mv "warmstart.dat" "RH80_PR5_10day_startfrom_4day_final.dat"
                
                elif [ -f "RH80_PR5_10day_startfrom_4day_final.dat" ] && [ ! -f "RH80_PR5_10day_startfrom_6day_final.dat" ]; then
                        mv "warmstart.dat" "RH80_PR5_10day_startfrom_6day_final.dat"
		
                elif [ -f "RH80_PR5_10day_startfrom_6day_final.dat" ] && [ ! -f "RH80_PR5_10day_startfrom_8day_final.dat" ]; then
                        mv "warmstart.dat" "RH80_PR5_10day_startfrom_8day_final.dat"
               		echo "All files have completed!!!"
			break	 
		else 
                        mv "warmstart.dat" "RH80_PR5_10day_startfrom_0day_final.dat"
		fi
		### run the next file ###
		julia Run_HS.jl 
	else
		echo 'The file does not exist.'
		break
	fi
done
