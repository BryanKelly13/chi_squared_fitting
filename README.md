To run this set of scripts

1st - Get a fresco output file from running fresco of desired reaction cross-section information,take the fort.16 file and move into
fresco output files directory. Next, run python3 getFresco_output.py <fort.16> and this will extract just the cross-section info
into a 2 column format. You also need your experimental cross-section values in a file w/ a 3 colum format - |labangle | cross-sec | err |
in the ang_dist folder.

2nd - In the main directory, run python3 merge.py <excitation energy>, where <excitation energy> is a search string to find all fresco
output files that match that energy, typically my naming scheme is finalnucleus_energy_transferconfiguration.txt 
An example would be for a 50Ti, 4172 keV state w/ l=1 transfer 2p3/2, the file name would be 50Ti_4172_2p32.txt, then the
getFresco_output.py will create 50Ti_4172_2p32.sorted

If you have different configurations you can try them separately, but do not have them both in the directory at once, b/c the current
working version will assume you have a mixed transfer (i.e. l=1+3) and try to create a file for this mixed transfer

3rd - Once you ran merge.py this will create an .inp file in the minimzing_folder directory, all you have left to do is in the main
directory run python3 chi_squared.py and this should output a plot with your DWBA curve scaled to data.
