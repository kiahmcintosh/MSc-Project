from mol_networking import Spectrum
import math


def test_normalisation():
    spectrum=Spectrum.Spectrum()
    spectrum.add_peak(1.0,605)
    spectrum.add_peak(2.0,234)
    spectrum.add_peak(3.0,574)
    spectrum.add_peak(4.0,348)
    spectrum.add_peak(5.0,943)
    spectrum.euclidean_scale()

    total=0
    for p in spectrum.peaks:
        total+=(p.scaled_intensity*p.scaled_intensity)
    assert math.isclose(total,1), "Incorrect normalisation"

