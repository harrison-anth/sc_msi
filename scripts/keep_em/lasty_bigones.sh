while read line; do if [[ ! -f /data3/hanthony/infer_cnv_temp/patient_$line/run.final.infercnv_obj ]]; then if [[ $line == "CRC2821" ]]; then Rscript infer_CRC2821.R /data3/hanthony/infer_cnv_temp/patient_CRC2821/ 110; else Rscript infer_patient_cnv.R $line /data3/hanthony/infer_cnv_temp/patient_$line/ 110; fi; fi; done < ../manifests/patient_ids5.txt

