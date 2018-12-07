#!/usr/bin/bash

# ./run_pipeline.sh TCGAcesc_adeno_squamous adeno squamous
# ./run_pipeline.sh TCGAthca_mutKRAS_mutBRAF mutKRAS mutBRAF

start_time=$(date -R)    
set -e

runDir="/mnt/etemp/marie/TAD_DE_pipeline_v2_TCGAdata"


TAD_DE_pipDir="/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom"
TAD_DE_script="./zzz_run_given_step_given_data_v2.sh"

function mvBack {
  echo "... go back to my folder"
  cd $runDir  
}
trap mvBack EXIT


echo_and_launch() {
   echo "> $1"
   $1
}


dataset="$1"
cond1="$2"
cond2="$3"



echo "*** START ***"
echo "... > dataset: $dataset"
echo "... > cond1: $cond1"
echo "... > cond2: $cond2"

datasetDir="/mnt/ed4/marie/other_datasets/$dataset"
settingF_outputFolder="/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput"
fpkmFile="$datasetDir/fpkmDT.Rdata"
rnaseqFile="$datasetDir/rnaseqDT_v2.Rdata"
cond1SampFile="$datasetDir/${cond1}_ID.Rdata"
cond2SampFile="$datasetDir/${cond2}_ID.Rdata"

#********************** HARD-CODED SETTINGS ********************************************

nPermut="10000"
ncpu="20"

minCpmRatio="20/888"

step1=1     # prepare setting file
step2=1     # run the pipeline


TAD_DE_pipSteps=( "0cleanInputTCGA" "1cleanInput" "2" "3" "5" "4" "6" "7" "8c" "9" "10" "11" "13cleanInput" "14f2" "170revision2EZH2" "170" )
#TAD_DE_pipSteps=( "0cleanInput" "1cleanInput" "2" "3" "5" )
#TAD_DE_pipSteps=( "0cleanInput" )
#TAD_DE_pipSteps=( "13cleanInput" )
#TAD_DE_pipSteps=( "170" )
#****************************************************************************************


new_setting_file="$settingF_outputFolder/run_settings_${dataset}.R"


if [[ "$step1" -eq 1 ]] ; then

	cat > ${new_setting_file} <<- EOM

    # > file written: `date -R` 

    # in this file, settings that are specific for a run on a dataset

    # gives path to output folder
    pipOutFold <- "OUTPUT_FOLDER/${dataset}"

    # full path (starting with /mnt/...)
    # following format expected for the input
    # colnames = samplesID
    # rownames = geneID
    # !!! geneID are expected not difficulted

    # *************************************************************************************************************************
    # ************************************ SETTINGS FOR 0_prepGeneData
    # *************************************************************************************************************************

    # UPDATE 07.12.2018: for RSEM data, the "analog" FPKM file is provided separately (built in prepData)
    rna_fpkmDT_file <- "$fpkmFile"

    rnaseqDT_file <- "$rnaseqFile"
    my_sep <- "\t"
    # input is Rdata or txt file ?
    # TRUE if the input is Rdata
    inRdata <- TRUE

    # can be ensemblID, entrezID, geneSymbol
    geneID_format <- "entrezID"
    stopifnot(geneID_format %in% c("ensemblID", "entrezID", "geneSymbol"))

    # are geneID rownames ? -> "rn" or numeric giving the column
    geneID_loc <- "rn"
    stopifnot(geneID_loc == "rn" | is.numeric(geneID_loc))

    removeDupGeneID <- TRUE

    # *************************************************************************************************************************
    # ************************************ SETTINGS FOR 1_runGeneDE
    # *************************************************************************************************************************

    # labels for conditions
    cond1 <- "$cond1"
    cond2 <- "$cond2"

    # path to sampleID for each condition - should be Rdata ( ! sample1 for cond1, sample2 for cond2 ! )
    sample1_file <- "$cond1SampFile"
    sample2_file <- "$cond2SampFile"

    minCpmRatio <- ${minCpmRatio} 

    inputDataType <- "RSEM"

    nCpu <- ${ncpu}
	
    # number of permutations
	nRandomPermut <- ${nPermut}
	step8_for_permutGenes <- TRUE
	step8_for_randomTADsFix <- FALSE
	step8_for_randomTADsGaussian <- FALSE
	step8_for_randomTADsShuffle <- FALSE
	step14_for_randomTADsShuffle <- FALSE

	EOM
	echo "WRITTEN: ${new_setting_file}"


echo "> END STEP1:" $(date -R)

fi # end if STEP1

###################################################################################################################################################
### STEP2: MOVE TO THE TAD DE PIPELINE DIRECTORY TO LAUNCH THE TAD DE PIPELINE
###################################################################################################################################################
if [[ "$step2" -eq 1 ]] ; then
	echo "> START STEP2:" $(date -R)
	echo_and_launch "cd $TAD_DE_pipDir"

	echo $TAD_DE_script ${new_setting_file} `echo ${TAD_DE_pipSteps[*]}`
	$TAD_DE_script ${new_setting_file} `echo ${TAD_DE_pipSteps[*]}`

	cd $runDir
	echo "> END STEP2:" $(date -R) 
fi


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

