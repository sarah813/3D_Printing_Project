# T3 EDS analysis

The PAD project contains two Axia ChemiSEM EDS map acquisitions at nominally the
same field of view: **T3-x1000-5kV** and **T3x1000-20kV**.

## Main interpretation

- Copper is the dominant metallic signal. The 20 kV spectrum provides the most
  defensible confirmation because it excites both Cu-L and the separated Cu-K
  lines near 8.05 and 8.90 keV.
- Carbon and oxygen are substantial low-energy signals. Carbon can include the
  specimen binder, mounting medium, carbon coating, or chamber contamination;
  oxygen is consistent with oxidized copper and/or oxygen-containing residue.
- In the maps, Cu is broadly distributed across the matrix, whereas O is
  concentrated around dark curved features and inclusions. This spatial
  anticorrelation is consistent with oxygen-rich non-matrix regions or surface
  oxidation, but the map alone cannot distinguish those mechanisms.
- The project also contains maps for N, Al, Si, P, S, Cl, Ca, and Ti. These
  should be treated as minor/trace candidates until their peaks and spatial
  distributions are checked against the spectrum and background.
- Cr, Nd, and La were selected automatically only in the 5 kV acquisition.
  Their low-energy L/M lines overlap strongly with O-K and Cu-L. They are not
  considered confirmed from this PAD alone and should be verified at higher
  voltage or with standards before reporting them as composition.
- The 20 kV acquisition has 44.4% detector dead time. It confirms the Cu-K
  peaks well, but precise quantitative comparison should use a lower count rate
  or vendor dead-time-corrected quantification.

## Quantification limitation

The PAD stores summed spectra and quantitative-map preview images, but it does
not contain a ready-made bulk weight-percent table. Preview intensity scales
are exported in `quant_map_channels.csv`; they are rendering scales, not bulk
composition values. Reliable wt%/at% should be exported from the vendor
software or recalculated with standards/background corrections.
