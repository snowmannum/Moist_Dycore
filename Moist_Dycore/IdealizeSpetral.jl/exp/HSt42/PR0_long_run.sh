#!bin/sh


start_day=0
final_day=20000
space_day=25

for i in `seq $start_day $space_day $final_day`
do
L=0
echo $i"day"
echo -n $L > Latent_heat.txt


if [ $i -eq 0 ]; then
	rm -rf HSt42_${L}
	mkdir HSt42_${L}
	echo -n $space_day > HSt42_${L}/day_interval.txt
	echo -n "None" > HSt42_${L}/firstday_file.txt
else
	echo -n $space_day > HSt42_${L}/day_interval.txt
	echo -n "warmstart_${L}.dat" > HSt42_${L}/firstday_file.txt
fi

julia Run_HS.jl
L=0 # To make sure that it wouldn't run the wrong L !!!
     # When there are many L run simultaneously, L might be mistaken.
     # For example, 
     # If you run the first file L = 0, it takes you about 1 hour to run julia Run_HS.jl, 
     # in this hour, if you run the second file, for example, L = 10, 
     # then after the first file finish julia Run_HS.jl,
     # it run the below code, it would take the output file to HSt42_10 (should've taken to HSt42_0) !!!
     # So to fix this, it should command second time to tell it the correct L !!!

if [ -f "HSt42_${L}/warmstart_${L}.dat" ] && [ $i -lt $final_day ]; then
	cp "HSt42_"${L}"/warmstart_${L}.dat" "warmstart_${L}.dat"
	mv "HSt42_"${L}"/all_L"${L}".dat" "HSt42_"${L}"/RH80_PR"$L"_"$final_day"day_startfrom_"$i"day_final.dat"
	echo 'warmstart file exists.'
elif [ -f "HSt42_${L}/warmstart_${L}.dat" ] && [ $i -eq $final_day ]; then
	mv "HSt42_"${L}"/all_L"${L}".dat" "HSt42_"${L}"/RH80_PR"$L"_"$final_day"day_startfrom_"$i"day_final.dat"
	mv "HSt42_"${L}"/warmstart_${L}.dat" "HSt42_"${L}"/HSt42_"${L}"RH80_PR"$L"_"$final_day"day_startfrom_"$i"day_final.dat"
	echo "All files have completed!!!"
fi
done
