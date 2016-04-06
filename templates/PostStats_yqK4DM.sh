#!/bin/bash

cd /home/arjen/out/
echo "Start poststats	" `date` "	" `uname -n` >> /home/arjen/out//logs/out.log

perl /opt/bamMetrics/bamMetrics.pl -bam /home/arjen/out//TESTY/mapping/TESTY_dedup.bam -bam /home/arjen/out//TESTZ/mapping/TESTZ_dedup.bam -bam /home/arjen/out//TESTX/mapping/TESTX_dedup.bam -output_dir /home/arjen/out//QCStats/ -run_name out -genome /data/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta -queue all.q -queue_threads 4 -queue_mem 16 -queue_time 256:0:0 -queue_project hartwig -picard_path /opt/picard-tools -debug -wgs -coverage_cap 250 
qalter -hold_jid bamMetrics_report_out,PostStats_yqK4DM PostStats_Check_74RUNB

