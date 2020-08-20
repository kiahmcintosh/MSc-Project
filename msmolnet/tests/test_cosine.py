"""Tests for similarity scoring methods
Use test_cosine() to run all tests
Throws an assertion error if a score or number of matched peaks is not correct
"""

from msmolnet import similarity
from msmolnet import Spectrum
import math

def same_spectrum():
    """Test two example spectra with same values
    """

    S1=Spectrum.Spectrum()
    S1.add_peak(50.7,234)
    S1.add_peak(54.6,585)
    S1.add_peak(60.7,773)
    S1.add_peak(65.6,387)
    S1.add_peak(87.7,546)
    S1.add_peak(104.6,598)
    S1.pep_mass=100
    S1.euclidean_scale()

    S2=Spectrum.Spectrum()
    S2.add_peak(50.7,234)
    S2.add_peak(54.6,585)
    S2.add_peak(60.7,773)
    S2.add_peak(65.6,387)
    S2.add_peak(87.7,546)
    S2.add_peak(104.6,598)
    S2.pep_mass=100
    S2.euclidean_scale()
   
    score,peaks=similarity.cosine_score_max(S1,S2)
    assert peaks==6, "Incorrect number of peaks matched with greedy method"
    assert math.isclose(score,1.0), "Incorrect score with greedy method"

    score,peaks=similarity.cosine_score_greedy(S1,S2)
    assert peaks==6, "Incorrect number of peaks matched with maximum weighted method"
    assert math.isclose(score,1.0), "Incorrect score with maximum weighted method"

def all_match():
    """Spectra where all the peaks match within the default tolerance
    """
    S1=Spectrum.Spectrum()
    S1.add_peak(50.7,234)
    S1.add_peak(54.6,585)
    S1.add_peak(60.7,773)
    S1.add_peak(65.6,387)
    S1.add_peak(87.7,546)
    S1.add_peak(104.6,598)
    S1.pep_mass=100
    S1.euclidean_scale()

    S2=Spectrum.Spectrum()
    S2.add_peak(50.5,234/2)
    S2.add_peak(54.8,585/2)
    S2.add_peak(61.0,773/2)
    S2.add_peak(65.4,387/2)
    S2.add_peak(88.0,546/2)
    S2.add_peak(104.3,598/2)
    S2.pep_mass=100
    S2.euclidean_scale()

    score,peaks=similarity.cosine_score_max(S1,S2)
    assert peaks==6, "Incorrect number of peaks matched with greedy method"
    assert math.isclose(score,1.0), "Incorrect score with greedy method"

    score,peaks=similarity.cosine_score_greedy(S1,S2)
    assert peaks==6, "Incorrect number of peaks matched with maximum weighted method"
    assert math.isclose(score,1.0), "Incorrect score with maximum weighted method"


def no_match():
    """Spectra where no peaks match
    """
    S1=Spectrum.Spectrum()
    S1.add_peak(50.7,234)
    S1.add_peak(54.6,585)
    S1.add_peak(60.7,773)
    S1.add_peak(65.6,387)
    S1.add_peak(87.7,546)
    S1.add_peak(104.6,598)
    S1.pep_mass=100
    S1.euclidean_scale()

    S2=Spectrum.Spectrum()
    S2.add_peak(50.2,234)
    S2.add_peak(53.8,585)
    S2.add_peak(61.3,773)
    S2.add_peak(66.2,387)
    S2.add_peak(88.1,546)
    S2.add_peak(103.9,598)
    S2.pep_mass=100
    S2.euclidean_scale()

    score,peaks=similarity.cosine_score_max(S1,S2)
    assert peaks==0, "Incorrect number of peaks matched with greedy method"
    assert score==0, "Incorrect score with greedy method"
    

    score,peaks=similarity.cosine_score_greedy(S1,S2)
    assert peaks==0, "Incorrect number of peaks matched with maximum weighted method"
    assert score==0, "Incorrect score with maximum weighted method"

def max_v_greedy():
    """Spectra which give a different score for greedy vs maximum-weighted
    """

    S1=Spectrum.Spectrum()
    S1.add_peak(50.4,16)
    S1.add_peak(50.7,36)
    S1.add_peak(74.8,25)
    S1.add_peak(96.2,23)
    S1.pep_mass=100
    S1.euclidean_scale()

    S2=Spectrum.Spectrum()
    S2.add_peak(50.6,49)
    S2.add_peak(50.9,25)
    S2.add_peak(74.6,9)
    S2.add_peak(102.4,17)
    S2.pep_mass=100
    S2.euclidean_scale()

    score,peaks=similarity.cosine_score_max(S1,S2)
    g_score,g_peaks=similarity.cosine_score_greedy(S1,S2)

    assert score>=g_score, "Maximum weighted method did not get higher score than greedy method"
    assert peaks>=g_peaks, "Maximum weighted method did not match more peaks than greedy method"

    assert peaks==3, "Incorrect number of peaks matched with greedy method"
    assert math.isclose(score,0.73), "Incorrect score with greedy method"

    assert g_peaks==2, "Incorrect number of peaks matched with maximum weighted method"
    assert math.isclose(g_score,0.57), "Incorrect score with maximum weighted method"

def test_cosine():
    """Run all cosine tests
    """
    same_spectrum()
    all_match()
    no_match()
    max_v_greedy()