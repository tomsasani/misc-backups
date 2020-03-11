set -euo pipefail
pwd
cd $TMPDIR

wget -O /usr/bin/goleft https://github.com/brentp/goleft/releases/download/v0.1.19/goleft_linux64
chmod +x /usr/bin/goleft

aws s3 sync s3://ql-tmp/crais .
aws s3 cp s3://ql-tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai .

goleft indexcov -d indexcov-out/ \
        --fai GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
            --sex "chrX,chrY" \
                *.crai

aws s3 sync --exclude "*.bed.gz" indexcov-out s3://ql-tmp/$name-indexcov/
