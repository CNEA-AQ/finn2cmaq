#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
#------------------------------------------------
#Input-Data:
start_date="2022-10-19"	#"%Y-%m-%d %H"
  end_date="2022-10-20"	#"%Y-%m-%d %H"

chemistry="GEOSchem" #GEOSchem # MOZ4 # SAPRC99

finn_data_directory="./finn_data/"
#------------------------------------------------
thisYear=$( date +"%Y" )

start_date_s=$(date -d "$start_date" +%s)
  end_date_s=$(date -d "$end_date  " +%s)

YYYY=$(date -d @$start_date_s +"%Y")   #año
 DDD=$(date -d @$start_date_s +"%j")   #dia juliano

day=$start_date_s
while [ $day -le $end_date_s ]
do
        YYYY=$(date -d @$day +"%Y")   #año
          MM=$(date -d @$day +"%m")   #mes
         DDD=$(date -d @$day +"%j")   #dia juliano
          HH=$(date -d @$day +"%H")   #hora
        echo "Dia: $YYYY-$DDD (mes: $MM)"


	#Get FINN data:
	 if [ ! -f $finn_data_directory/GLOB_${chemistry}_$YYYY$DDD.txt ]
	 then
	 	if [ $YYYY -lt $(($thisYear-2)) ]	#los años viejos los tienen en carpetas.
	 	then
	 	       wget "https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/$YYYY/GLOB_${chemistry}_$YYYY$DDD.txt.gz" -P finn_data/
	 	else
	 	       wget "https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/GLOB_${chemistry}_$YYYY$DDD.txt.gz" -P finn_data/
	 	fi;
		
		gzip -d finn_data/GLOB_${chemistry}_$YYYY$DDD.txt.gz
	fi;
	finnFile=finn_data/GLOB_${chemistry}_$YYYY$DDD.txt


        day=$(($day + 86400 )) # Incrementar day un dia
done

