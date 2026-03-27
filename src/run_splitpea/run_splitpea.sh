# install splitpea through https://github.com/ylaboratory/splitpea

## 1
python src/combine_spliced_exons.py test-data

## 2
Rscript src/delta_psi.R -o psis -s test-data/spliced_exons_gtex_pancreas_test_combined_mean.txt -b test-data/spliced_exons_gtex_pancreas_test.txt -t test-data/spliced_exons_tcga_paad_test.txt
#Rscript src/delta_psi.R -o psis -s test-data/gtex_combined_mean.txt -b test-data/gtex.txt -t test-data/tcga.txt

## 3
python src/get_background_ppi.py

## 4
bash src/run_splitpea_batch.sh -i psis -o /results/splitpea_outputs -p 75

## 5
python src/get_consensus_network.py /results/splitpea_outputs