# chrono-python
This is the github repository for wave-fitting for chronobiology.

## Sine wave fitting 
Python requirements (run the following):

`conda install numpy scipy pandas matplotlib jupyter`

## How is the data for fitting supposed to be formatted?
For every gene...

Create an excel sheet with .csv file format. For a block, 4 columns (`Cq_c SEM_c Cq_e SEM_e`), and 6 rows (for 6 zeitgeber time readings).

If multiple blocks are present, place them one below the other, with a row of (`1e-5 1e-5 1e-5 1e-5`) separating the blocks.

