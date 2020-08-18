from msmolnet.use_matchms import convert_matches as convert

import numpy as np
import matchms
from matchms import Spectrum
from matchms.similarity import CosineGreedy



def test_convert_matches():
    spectrum = Spectrum(mz=np.array([100, 150, 200.]),
                intensities=np.array([0.7, 0.2, 0.1]),
                metadata={'SCANS': '1'})

    new_spectrum=convert.convert_spectrum(spectrum)

    assert new_spectrum.peaks[0].mass==100 and new_spectrum.peaks[0].intensity==0.7, "incorrect peak"
    assert new_spectrum.feature_id=='1', "incorrect ID"