# slim tree
*You have to provide a configuration file containing input file paths for every dataset (and additional info) + the short dataset name to analyze (defined in the configuration file)*

python3 python/slim_hnl_tree.py cfg/*configuration_file_name.json* *short_dataset_name*

*e.g.*
python3 python/slim_hnl_tree.py cfg/hnl_tree_input_fromCrab.json QCD_Pt-20to30

#analyze tree

*You have to provide the same info as above*

python3 python/hnl_tree_analyzer.py cfg/*configuration_file_name.json* *short_dataset_name*

*e.g.*
python3 python/hnl_tree_analyzer.py cfg/hnl_tree_input_fromCrab.json QCD_Pt-20to30

#make data VS. mc plots

python3 python/make_dataVSmc_comparison.py dataVSmc_input.json

