#!/bin/bash -l
#$ -l h_vmem=250M

###
# Cluster array script for running multiple IKs/IKr ratios. Array should have
# have range of 0-40.
###

# This script is meant to be run as an array job. SGE_TASK_ID should exist if
# this is true.
# Using old-school POSIX style check since cluster uses Bash 4.1
if [ -n "${SGE_TASK_ID+1}" ]; then
    echo "Array job detected: SGE_TASK_ID = $SGE_TASK_ID"
else
    echo "ERROR: SGE_TASK_ID does not exist, check if job was run as array"
    exit 1
fi

# Conductance values are hard coded.
IKsArray=(0 0.0228322847894584 0.0413800525225716 0.0569770769295159 \
0.0702730282089285 0.0822805377584615 0.0924146710078945 0.101641057233047 \
0.110468256203516 0.117988205448277 0.124863475239707 0.131757060193437 \
0.137468142949850 0.143161380047815 0.148351379656497 0.151807193953271 \
0.156177064948843 0.160536709450943 0.164667901774919 0.167886745111491 \
0.172360954499634 0.175605820186594 0.178274517665887 0.180911898084707 \
0.184327311571487 0.186255036846788 0.189743769086285 0.191392633213498 \
0.194230425893375 0.197046372546094 0.198447404077553 0.200844299268327 \
0.202975393916574 0.204502991266026 0.206645695508080 0.208103819538936 \
0.210352084984118 0.211561226902314 0.213690412958983 0.214666251167852 \
0.215987849297908)

IKrArray=(0.0512653261862278 0.0456645695789169 0.0413800525225716 \
0.0379847179530106 0.0351365141044643 0.0329122151033846 0.0308048903359648 \
0.0290403020665849 0.0276170640508789 0.0262196012107283 0.0249726950479414 \
0.0239558291260794 0.0229113571583083 0.0220248276996638 0.0211930542366425 \
0.0202409591937694 0.0195221331186054 0.0188866717001110 0.0182964335305465 \
0.0176722889591043 0.0172360954499634 0.0167243638272946 0.0162067743332624 \
0.0157314693986702 0.0153606092976239 0.0149004029477430 0.0145956745450988 \
0.0141772320898888 0.0138736018495268 0.0135894050031789 0.0132298269385035 \
0.0129576967269889 0.0126859621197859 0.0123941206827895 0.0121556291475341 \
0.0118916468307963 0.0116862269435621 0.0114357419947197 0.0112468638399465 \
0.0110085257009155 0.0107993924648954)

IDX=$SGE_TASK_ID

echo "Starting Restitution Protocol Simulation: Run $IDX"
echo "IKs: ${IKsArray[$IDX]}"
echo "IKr: ${IKrArray[$IDX]}"

# Load GCC compiler
slchoose gcc 4.7.4 gcc4_64

# Move into node's temp directory
cd $TMPDIR

# Copy repository and compile
rsync -av ${HOME}/Active_Development/LivRudy2009 .
cd LivRudy2009/exe/IKs_IKr_Restitution
make clean
make
echo "Finished copying repository and compiling"

# Script variables
SIM_EXEC=${TMPDIR}/LivRudy2009/exe/IKs_IKr_Restitution/IKs_IKr_Restitution
DATA_DIR=${HOME}/Data/2017/LivRudy2009
DATA_FOLDER=IKs_IKr_Restitution
DATA_FILE=IKs_IKr_Data_$IDX.dat
LOG_FILE=IKs_IKr_Data_$IDX.log

chmod a+x $SIM_EXEC

# Run simulation
echo $SIM_EXEC ${IKsArray[$IDX]} ${IKrArray[$IDX]} $DATA_FILE
$SIM_EXEC ${IKsArray[$IDX]} ${IKrArray[$IDX]} $DATA_FILE 1>> $LOG_FILE
echo "Simulation complete"

echo "Copying files to home directory"
# Rsync data into specified data directory
# First make sure directory exits, if not create it
if [ ! -d "$DATA_DIR/$DATA_FOLDER" ]; then
    mkdir -p $DATA_DIR/$DATA_FOLDER
fi
# Rsync data and log file
rsync -av $DATA_FILE $DATA_DIR/$DATA_FOLDER
rsync -av $LOG_FILE $DATA_DIR/$DATA_FOLDER

echo "Completed restitution protocol simulation"

# Successful completion
exit 0
