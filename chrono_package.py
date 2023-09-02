import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize

class CHRONO():
    
    def __init__(self, address=None, NN=6, block=0, name=''):
        
        self.NN       = NN
        self.block    = block
        self.address  = address
        self.name     = name

        self.tt  = np.linspace(2,22,self.NN)
        self.tt2 = np.linspace(2,22,50)
        
        print('\033[1m', 'Analyzing {} gene...\n'.format(self.name))
#         print(self.address)


    def get_measurements(self):
        print('\033[1m', 'Reading the csv file...\n')
        '''Provide address of your excel sheet containing the block data, number of readings in a block, and the block number'''

        sheet = pd.read_csv(self.address)

        li,ri = (self.NN+1)*(self.block), (self.NN+1)*(self.block)+self.NN

        self.cqc = sheet['Cq_c'][li:ri]
        self.erc = sheet['SEM_c'][li:ri]

        self.cqe = sheet['Cq_e'][li:ri]
        self.ere = sheet['SEM_e'][li:ri]


        return self.cqc, self.erc, self.cqe, self.ere


    def get_sin_params(self, data, error, limits=[20,28], print_guess_and_fitted=True):
#         print('\033[1m', 'Obtain the parameters for sine wave fitting...')
        
        '''Provide DATA and ERROR, limits of periodicity (default is set as [20,28] hours)'''

        data  = np.array(data)
        ff    = np.fft.fftfreq(len(self.tt), (self.tt[1]-self.tt[0]))   
        Fdata = abs(np.fft.fft(data))

        guess_freq   = abs(ff[np.argmax(Fdata[1:])+1])   

        if 1/guess_freq<limits[0]:
            guess_freq = 1/limits[0]

        if 1/guess_freq>limits[1]:
            guess_freq = 1/limits[1]

        guess_amp    = np.std(data) * 2.**0.5
        guess_offset = np.mean(data)
        guess        = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

        def sinfunc(t, A, w, p, c):  
#             print('sinfunc = ', A * np.sin(w*t + p) + c)
            return A * np.sin(w*t + p) + c

        bounds = ([-1e3, 2*np.pi*(1/limits[1]), -1e3, -1e3], [1e3, 2*np.pi*(1/limits[0]), 1e3, 1e3])
        
#         print(guess, bounds, error)

        self.popt, self.pcov = optimize.curve_fit(sinfunc, self.tt, data, p0=guess, sigma=error, bounds=bounds)
        A, w, p, c           = self.popt
        f                    = w/(2*np.pi)
        self.sigmas          = np.sqrt(np.diag(self.pcov))
        
        if print_guess_and_fitted:
            print('Guessed period = {:.1f} hrs.'.format(1/guess_freq))
            print('Fitted period  = {:.1f} hrs.'.format(abs(1/f)))
            print('\n')
            
        return self.popt, self.sigmas
    
    
    
    
        
    def fit_sinwave(self, data, error):
#         print('\033[1m', 'Obtain the sine wave...')
        
        "Provide DATA and ERROR to obtain the fitted sine wave and the error in fit."
        
        [A, w, p, c], [Ae, we, pe, ce]  = self.get_sin_params(data, error)
        
#         self.sinusoid    = lambda t: A * np.sin(w*t + p) + c

        self.sinusoid    = A * np.sin(w*self.tt2 + p) + c
    
        angle = w*self.tt2 + p 
        df_dA =                np.sin(angle) * Ae
        df_dw = A * self.tt2 * np.cos(angle) * we
        df_dp =            A * np.cos(angle) * pe
        df_dc =                            1 * ce

        self.df = np.sqrt(df_dA**2 + df_dw**2 + df_dp**2 + df_dc**2)
          
        return self.sinusoid, self.df




    def plot_curves(self, DATA, plot_range=[-0.06,0.06], title='', save_address='', which='both'):
        
        '''Provide DATA (and ERROR), plot range for y-axis, plot title (optional)'''

        print('\033[1m', '\n Plot the fitted sine wave...\n')
        
        plt.figure(figsize=(5,4), dpi=120)
        plt.title(title)
        
        if which=='control':
            ct, error_ct = self.fit_sinwave(DATA[0], DATA[1])
            plt.scatter(self.tt, self.cqc, marker='x', color='r', s=80)
            plt.errorbar(self.tt, self.cqc, self.erc, fmt='none', color='r')
            plt.plot(self.tt2, ct, "r--", label="control", linewidth=2)
            plt.fill_between(self.tt2, ct+error_ct, ct-error_ct, color='red', alpha=0.4)

        elif which=='experimental':    
            ex, error_ex = self.fit_sinwave(DATA[2], DATA[3])
            plt.scatter(self.tt, self.cqe, marker='d', color='b', s=80)
            plt.errorbar(self.tt, self.cqe, self.ere, fmt='none', color='b')
            plt.plot(self.tt2, ex, "b:", label="experimental", linewidth=2)
            plt.fill_between(self.tt2, ex+error_ex, ex-error_ex, color='blue', alpha=0.4)

        else:
            ct, error_ct = self.fit_sinwave(DATA[0], DATA[1])
            ex, error_ex = self.fit_sinwave(DATA[2], DATA[3])
            plt.scatter(self.tt, self.cqc, marker='x', color='r', s=80)
            plt.errorbar(self.tt, self.cqc, self.erc, fmt='none', color='r')
            plt.plot(self.tt2, ct, "r--", label="control", linewidth=2)
            plt.fill_between(self.tt2, ct+error_ct, ct-error_ct, color='red', alpha=0.4)
            
            plt.scatter(self.tt, self.cqe, marker='d', color='b', s=80)
            plt.errorbar(self.tt, self.cqe, self.ere, fmt='none', color='b')
            plt.plot(self.tt2, ex, "b:", label="experimental", linewidth=2)
            plt.fill_between(self.tt2, ex+error_ex, ex-error_ex, color='blue', alpha=0.4)

            
        plt.xticks(self.tt)
        plt.xlabel('time, hrs.')
        plt.grid(axis='both')
        plt.legend()
        
        if plot_range=='smart':
            dat = np.array(DATA)
            plot_range = [2*dat.min(), 2*dat.max()]
        
        plt.ylim(plot_range[0], plot_range[1])
        plt.xlim(0,24)
        
        print('\033[3m', 'Saving the plot in file named {}_{:02d}.png'.format(self.name, self.block))
        plt.savefig('{}_{:02d}.pdf'.format(self.name, self.block))
        plt.show()





print('\033[1m', 'IF YOU \U0001FAF5 ARE RAISED RIGHT, YOU WILL GIVE PRASAD AUTHORSHIP \N{flexed biceps}. YOU ARE WELCOME. \U0001F60E')




















