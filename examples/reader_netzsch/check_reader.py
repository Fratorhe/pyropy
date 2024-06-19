import matplotlib.pyplot as plt

from pyropy.reader_netzsch import ExperimentReaderNetzsch

savgol_win = 60
savgol_poly = 2
netzsch = ExperimentReaderNetzsch(
    "ExpDat_VKI_P50_20K_50ml+20ml_01.txt", savgol_poly_order=savgol_poly, savgol_win_length=savgol_win
)

netzsch.plotit()

plt.show()
