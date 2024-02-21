#!/bin/bash
myArray=(drugpipeline.sh 1a_multivariable_MR.R 2a_observeresults.R 2b_plot_all_figures.R 2c_write_tables.R 2d_writetext.R)
for str in ${myArray[@]}; do
chmod u+x ./$str
done
echo 'Initializing drugpipeline.sh' && ./drugpipeline.sh && echo 'Initializing 1a_multivariable_MR.R' && ./1a_multivariable_MR.R && echo 'Initializing 2a_observeresults.R' && ./2a_observeresults.R && echo 'Initializing 2b_plot_all_figures.R' && ./2b_plot_all_figures.R && echo 'Initializing 2c_write_tables.R' && ./2c_write_tables.R && echo 'Initializing 2d_writetext.R' && ./2d_writetext.R && echo 'The master script finished without errors'
