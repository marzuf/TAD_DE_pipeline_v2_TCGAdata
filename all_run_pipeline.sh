#!/bin/bash

# ./run_pipeline.sh "TCGAcesc_adeno_squamous adeno squamous

start_time=$(date -R)    
set -e

all_args=(      
"TCGAbrca_lum_bas lum bas"
)

#all_args=(      
#"TCGAcoad_msi_mss msi mss"
#)
#all_args=(
#"TCGAcesc_adeno_squam adeno squam"
#"TCGAhnsc_HPVneg_HPVpos HPVneg HPVpos"
#"TCGAlgg_IDHwt_IDHmutnc IDHwt IDHmutnc"
#"TCGAsarc_ddlps_lms ddlps lms"
#"TCGAsarc_ddlps_mfs ddlps mfs"
#"TCGAsarc_lms_mfs lms mfs"
#"TCGAstad_EBVpos_gs EBVpos gs"
#"TCGAstad_msi_gs msi gs"
#"TCGAstad_norm_gs norm gs"
#"TCGAtgct_sem_nonsem sem nonsem"
#"TCGAucec_msi_cnl msi cnl"
#"TCGAacc_wt_mutCTNNB1 wt mutCTNNB1"
#"TCGAlihc_wt_mutCTNNB1 wt mutCTNNB1"
#"TCGAluad_wt_mutKRAS wt mutKRAS"
#"TCGApaad_wt_mutKRAS wt mutKRAS"
##"TCGAlaml_wt_mutFLT3 wt mutFLT3"
#"TCGAskcm_wt_mutBRAF wt mutBRAF"
#"TCGAskcm_wt_mutCTNNB1 wt mutCTNNB1"
#"TCGAthca_wt_mutBRAF wt mutBRAF"
#"TCGAblca_norm_blca norm blca"
#"TCGAkich_norm_kich norm kich"
#"TCGAlusc_norm_lusc norm lusc"
#"TCGAluad_mutKRAS_mutEGFR mutKRAS mutEGFR"
#"TCGAthca_mut.RAS_mutBRAF mut.RAS mutBRAF"
#)
# N=24

#all_args=(
#"TCGAlaml_wt_mutFLT3 wt mutFLT3"
#)

#all_args=(
#"TCGAgbm_classical_neural classical neural"
#"TCGAgbm_classical_proneural classical proneural"
#"TCGAgbm_classical_mesenchymal classical mesenchymal"
#"TCGAluad_nonsmoker_smoker nonsmoker smoker"
#)

#all_args=(
#"TCGAskcm_lowInf_highInf lowInf highInf"
#)

for arg in "${all_args[@]}"; do

    echo ./run_pipeline.sh $arg
    ./run_pipeline.sh $arg

done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

