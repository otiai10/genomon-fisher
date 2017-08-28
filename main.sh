
echo "[genomon-fisher][0] Checking input parameters (mostly provided as env variables)"
if [ ! -f "/var/refs/${REFERENCE}" ]; then
    echo "[genomon-fisher][ERROR] REFERENCE file not found: specified with '${REFERENCE}', searched at /var/refs"
    exit 1
fi
if [ ! -f "/var/data/${INPUT01}" ]; then
    echo "[genomon-fisher][ERROR] INPUT01 file not found: specified with '${INPUT01}', searched at /var/data"
    exit 1
fi
if [ ! -f "/var/data/${INPUT02}" ]; then
    echo "[genomon-fisher][ERROR] INPUT02 file not found: specified with '${INPUT02}', searched at /var/data"
    exit 1
fi

RESULT_PREFIX=result.${INPUT01}_x_${INPUT02}

echo "[genomon-fisher][1] >>> Decompress files if having \".tar.bz2\" extension: '${INPUT01}', '${INPUT02}'"
if [[ ${INPUT01} == *.tar.bz2 ]]; then
  tar xvjf /var/data/${INPUT01} --directory /var/data
  INPUT01=`echo ${INPUT01} | sed "s/\.tar.bz2$//"`
fi
if [[ ${INPUT02} == *.tar.bz2 ]]; then
  tar xvjf /var/data/${INPUT02} --directory /var/data
  INPUT02=`echo ${INPUT02} | sed "s/\.tar.bz2$//"`
fi

echo "[genomon-fisher][2] >>> align provided read samples: '${INPUT01}', '${INPUT02}'"
/bin/bwa mem /var/refs/${REFERENCE} /var/data/${INPUT01} /var/data/${INPUT02} > /var/data/${RESULT_PREFIX}.sam

echo "[genomon-fisher][3] >>> bamify result sam file: '${RESULT_PREFIX}.sam'"
/bin/samtools view -S -b -h /var/data/${RESULT_PREFIX}.sam > /var/data/${RESULT_PREFIX}.bam

echo "[genomon-fisher][4] >>> sort bam: '${RESULT_PREFIX}.bam'"
/bin/samtools sort /var/data/${RESULT_PREFIX}.bam > /var/data/${RESULT_PREFIX}.sorted.bam

echo "[genomon-fisher][5] >>> execute Fisher's test: '${RESULT_PREFIX}.sorted.bam'"
/usr/bin/fisher single \
	-o /var/data/${RESULT_PREFIX}.fisher.txt \
	--ref_fa /var/refs/$REFERENCE \
	-1 /var/data/${RESULT_PREFIX}.sorted.bam \
	--samtools_path /bin/samtools \
	--min_depth 8 --base_quality 15 --min_variant_read 4 --min_allele_freq 0.02 --post_10_q 0.02 \
	--samtools_params "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"

echo "[genomon-fisher][6] >>> Everything completed, bye."
