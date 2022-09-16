#1/bin/bash
# These example commands will run the targeted NOMe pipeline 
# from start to finish. Just modify the four fields below accordingly.
# Run within an interactive session with 8 cores:
# ish -l h_vmem=4G -pe smp 8 -binding linear:8

data="/seq/epiprod02/Battaglia/NanoNOMe/190522_GM12878/"
name="190522"
locus="190522_loci_tmp.bed"
dir="/seq/epiprod02/kdong/targetedNOMe_2022/pipeline"

cd $data

${dir}/1_primary_analysis.sh --index --align --filter --meth -n ${name} -L ${locus} \
	-f fast5/ > primary_analysis.log 2>&1

${dir}/2_visualize.sh --avg --smoo --bsm -n ${name} > \
	visualize.log 2>&1

${dir}/3_r_analysis.sh -n ${name} -L ${locus} > \
	r_analysis.log 2>&1

${dir}/4_haplotype.sh --var --bsm --avg --smoo --bw -n ${name} -L ${locus} > \
	haplotype.log 2>&1