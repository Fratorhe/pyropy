import matplotlib.pyplot as plt

from pyropy.reader_netzsch import ExperimentReaderNetzsch

savgol_win = 60
savgol_poly = 2
netzsch = ExperimentReaderNetzsch(
    "ExpDat_VKI_P50_20K_50ml+20ml_01.txt", savgol_poly_order=savgol_poly, savgol_win_length=savgol_win
)

netzsch.plotit()

plt.show()


# since the original Protocol for file reader only wants the filename and the folder, we can do a partial application of the filter

from functools import partial

# Create a partial function that always sets savgol_poly_order=savgol_poly and savgol_win_length=savgol_win
ExperimentReaderNetzschSavgol = partial(
    ExperimentReaderNetzsch, savgol_poly_order=savgol_poly, savgol_win_length=savgol_win
)

netzsch2 = ExperimentReaderNetzsch(
    "ExpDat_VKI_P50_20K_50ml+20ml_01.txt", savgol_poly_order=savgol_poly, savgol_win_length=savgol_win
)

netzsch2.plotit()

plt.show()
