# chrono-python
This is the github repository for wave-fitting for chronobiology.

## Sine wave fitting 
Python requirements (run the following):

`conda install numpy scipy pandas matplotlib jupyter`

## How is the data for fitting supposed to be formatted?
For every gene...

Create an excel sheet with .csv file format. A block must contain 4 columns (`Cq_c SEM_c Cq_e SEM_e`) and N rows (for N zeitgeber-time readings).

If multiple blocks are present, place them one below the other, with a row of (`1e-5 1e-5 1e-5 1e-5`) separating the blocks.

All you need to provide for the code to run is the name of the gene, address of the `.csv` file, block number, and number of zeitgeber-time readings.

## Functionalities in the code
- `get_measurements`: gives the `Cq_c, SEM_c, Cq_e, SEM_e`
- `get_sin_params`  : gives the parameters `[A,w,p,c]` found for fitting the sine wave `A * np.sin(w*t + p) + c`, and their respective errors `[Ae, we, pe, ce]`.
- `fit_sinwave`     : gives the sinewave and the error on sinewave, on a smooth array.
- `plot_curves`     : plots the data and their standard deviation, and the fitted sinewave and the error in fitting.

A demo code is provided in the notebook file named demo_chrono.ipynb

If you find this useful, please consider citing it in your work : `https://github.com/captain-shrugging-pants/chronobiology-python.git`


