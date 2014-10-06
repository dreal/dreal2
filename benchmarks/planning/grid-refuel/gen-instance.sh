#!/bin/bash

DIMENSION=$1

ADJACENCY="//set<pair<int, int>*> adjacent = {"


outputInit(){
    echo "init:"
    echo "@1 (and (fuel_level = 0) (distance_traveled = 0) (velocity = 0) (acceleration = 0) (clock = 0));"
}

outputGoal(){
    END=`expr ${DIMENSION} \* 4`
    END=`expr $END + 1`
    echo "goal:"
    echo "@${END} true;"
}


outputLocationMode(){
    MODE=$1
    REFUEL_MODE=`expr $MODE + 1`
    ACCEL_MODE=`expr $MODE + 2`
    DECEL_MODE=`expr $MODE + 3`
    MAX_MODE=`expr $DIMENSION \* 4`
    MAX_MODE=`expr $MAX_MODE + 1`
    
    echo "{ mode ${MODE};"

    echo "  flow:"
    echo "	d/dt[fuel_level] = 0;"
    echo "	d/dt[distance_traveled] = 0;"
    echo "	d/dt[velocity] = acceleration;"
    echo "	d/dt[acceleration] = 0;"
    echo "	d/dt[clock] = 1;"
    echo "  jump:"
    if [ $MODE -lt $MAX_MODE ] ; then
    echo "	//refuel"
    echo "	(and (clock >= epsilon) (velocity = 0)) ==> @"$REFUEL_MODE" (clock' = 0);"
    echo "	//drive"
    echo "	(and (clock >= epsilon)) ==> @"$ACCEL_MODE" (and  (clock' = 0));"
    else
	echo ""
#    echo "	true ==> @"$MODE" true;"
    fi
    echo "}"
}

outputRefuelMode(){
    MODE=$1
    REFUEL_MODE=`expr $MODE + 1`
    echo "{ mode ${REFUEL_MODE};"
    echo "  invt:"
    echo "      (clock <= 5);"
    echo "  flow:"
    echo "	d/dt[fuel_level] = station_"$MODE"_rate;"
    echo "	d/dt[distance_traveled] = 0;"
    echo "	d/dt[velocity] = 0;"
    echo "	d/dt[acceleration] = 0;"
    echo "	d/dt[clock] = 1;"
    echo "  jump:"
    echo "	//refuel"
    echo "	(clock >= epsilon) ==> @"$MODE" (clock' = 0);"
    echo "}"
}


outputAccelMode(){
    MODE=$1
    ACCEL_MODE=`expr $MODE + 2`
    DECEL_MODE=`expr $MODE + 3`
    NEXT_MODE=`expr $MODE + 4`
    echo "{ mode ${ACCEL_MODE};"

    echo "  flow:"
    echo "	d/dt[fuel_level] = -(0.1 * (velocity^2));"
    echo "	d/dt[distance_traveled] = velocity;"
    echo "	d/dt[velocity] = acceleration;"
    echo "	d/dt[acceleration] = 0;"
    echo "	d/dt[clock] = 1;"
    echo "  jump:"
    echo "	//drive"
    echo "	(and (clock >= epsilon)  (distance_traveled = distance_"$MODE"_"$NEXT_MODE")) ==> @"$NEXT_MODE" (and (distance_traveled' = 0) (clock' = 0));"
    echo "	//accel"
    echo "	(clock >= epsilon) ==> @"$ACCEL_MODE" (and (acceleration' = (acceleration + 1)) (clock' = 0));"
    echo "	//decel"
    echo "	(clock >= epsilon) ==> @"$DECEL_MODE" (and (acceleration' = (acceleration - 1)) (clock' = 0));"
    echo "}"
}

outputDecelMode(){
    MODE=$1
    ACCEL_MODE=`expr $MODE + 2`
    DECEL_MODE=`expr $MODE + 3`
    NEXT_MODE=`expr $MODE + 4`
    echo "{ mode ${DECEL_MODE};"

    echo "  flow:"
    echo "	d/dt[fuel_level] = -(0.1 * (velocity^2));"
    echo "	d/dt[distance_traveled] = velocity;"
    echo "	d/dt[velocity] = acceleration;"
    echo "	d/dt[acceleration] = 0;"
    echo "	d/dt[clock] = 1;"
    echo "  jump:"
    echo "	//drive"
    echo "	(and (clock >= epsilon)  (distance_traveled = distance_"$MODE"_"$NEXT_MODE")) ==> @"$NEXT_MODE" (and (distance_traveled' = 0) (clock' = 0));"
    echo "	//accel"
    echo "	(clock >= epsilon) ==> @"$ACCEL_MODE" (and (acceleration' = (acceleration + 1)) (clock' = 0));"
    echo "	//decel"
    echo "	(clock >= epsilon) ==> @"$DECEL_MODE" (and (acceleration' = (acceleration - 1)) (clock' = 0));"
    echo "}"
}


for((i=0;i<${DIMENSION};i++)); do {
	THIS=`expr $i \* 4`
	THIS=`expr $THIS + 1`
	NEXT=`expr $i + 1`
	NEXT=`expr $NEXT \* 4`
	NEXT=`expr $NEXT + 1`
	echo "#define distance_"$THIS"_"$NEXT" 10"
	echo "#define station_"$THIS"_rate 2"
} ; done
echo "#define epsilon 0.0001"

echo "[0, 100] fuel_level;"
echo "[0, 100] distance_traveled;"
echo "[0, 100] velocity;"
echo "[-5, 5] acceleration;"
echo "[0, 100] clock;"
echo "[0, 100] time;"

for((i=0;i<${DIMENSION};i++)); do {
	GROUP=`expr 4 \* $i`
	GROUP=`expr $GROUP + 1`
	outputLocationMode $GROUP
	outputRefuelMode $GROUP
	outputAccelMode $GROUP
	outputDecelMode $GROUP 
}; done
GROUP=`expr 4 \* $DIMENSION`
outputLocationMode `expr $GROUP + 1`


outputInit
outputGoal

