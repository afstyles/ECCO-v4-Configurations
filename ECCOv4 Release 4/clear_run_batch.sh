#!/bin/bash

# REMEMBER - No more than 16 jobs submitted on ARCHER at any time

# Set to none if you don't want to use any perturbations

declare -a PERTLABELS=("none"
                       # Qnet perturbations >>>>>>>>>>>>>>>>>>>>>
                       #"ArcticHeatBand_RAPID_upper_pulse_1997-03-01_1997-03-31_minus"
                       #BottomDrag perturbations
                       #"1e-4nonlin_plus"
                       #"1e-4nonlin_minus"
                     )

declare -a TMASKS=("DECend_1day_5year_mask"
                   "JUNend_1day_5year_mask"
                   #"MARend_1day_5year_mask"
                   #"SEPend_1day_5year_mask"
                  )

declare -a CMASKDIRS=(  #OSNAP (upper)
                        "OSNAP_upperExp"
                        "OSNAP_upperExp"
                        "OSNAP_upperExp"
                        "OSNAP_upperExp"
                        #OSNAP (lower)
                        #"OSNAP_lowerExp"
                        #"OSNAP_lowerExp"
                        #"OSNAP_lowerExp"
                        #"OSNAP_lowerExp"
                        # Global (global xy area weighted)
                        # "rA"
                      )

declare -a CMASKS=(   #OSNAP (upper)
                      "OSNAP_upperExp_East_100_1000_mask"
                      "OSNAP_upperExp_GRNE_100_1000_mask"
                      "OSNAP_upperExp_GRNW_100_1000_mask"
                      "OSNAP_upperExp_West_100_1000_mask"
                      #OSNAP (lower)
                      #"OSNAP_lowerExp_East_1000_3000_mask"
                      #"OSNAP_lowerExp_GRNE_1000_3000_mask"
                      #"OSNAP_lowerExp_GRNW_1000_3000_mask"
                      #"OSNAP_lowerExp_West_1000_3000_mask"
                      # Global (xy area weighted)
                      #"rA_mask"
                  )

for k in "${!PERTLABELS[@]}"; do
   for i in "${!CMASKS[@]}"; do
      for j in "${!TMASKS[@]}"; do

         echo "i = $i"
         echo "j = $j"
         echo "k = $k"

         TMASK=${TMASKS[j]}
         CMASK=${CMASKS[i]}
         CMASKDIR=${CMASKDIRS[i]}
         PERTLABEL=${PERTLABELS[k]}

         echo "TMASK = ${TMASK}"
         echo "CMASK = ${CMASK}"
         echo "CMASKDIR = ${CMASKDIR}"
         echo "PERTLABEL = ${PERTLABEL}"

         if [ "${PERTLABEL}" = none ] ; then
            export LABEL=${CMASK}.${TMASK}
         else
            export LABEL=${PERTLABEL}.${CMASK}.${TMASK}
         fi

         cp -v clear_run.template.sh clear_run.${LABEL}.sh
         sed -i "s/___TMASK___/${TMASK}/g" clear_run.${LABEL}.sh
         sed -i "s/___CMASK___/${CMASK}/g" clear_run.${LABEL}.sh
         sed -i "s/___CMASKDIR___/${CMASKDIR}/g" clear_run.${LABEL}.sh
         sed -i "s/___PERTLABEL___/${PERTLABEL}/g" clear_run.${LABEL}.sh
         sed -i "s/___LABEL___/${LABEL}/g" clear_run.${LABEL}.sh

         sh clear_run.${LABEL}.sh

         #Submit shell script
         sh runCommand.${LABEL}.sh
      done
   done
done
