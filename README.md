# Comment on "Inner Core Rotation Captured by Earthquake Doublets and Twin Stations" by Yang and Song

**Dongdong Tian** and **Lianxing Wen**

## Abstract

[Yang & Song (2022)](https://doi.org/10.1029/2022GL098393)
first claimed existence of Earth's inner core differential rotation
based on the waveform similarity of two neighboring stations AAK and KZA across an
earthquake doublet and then postulated a local velocity gradient at the top of the inner core
based on the difference of PKiKP-PKIKP differential times between the stations and
inferred inner core differential rotation rate. In this comment, we collectively analyze
the seismic data in the region and add the data of another nearby station HORS into analysis.
HORS and KZA, located in an opposite direction away from AAK, consistently exhibit
high waveform similarity. Collective analysis of seismic data demonstrates the invalidity
of both their logic of claiming existence of inner core differential rotation and
their postulation of "a local inner core gradient" to infer differential rotation.
Localized and episodic inner core surface change provides a physically consistent
explanation to the seismic data.

## Content of this repository

This repository contains the data and scripts used in Tian & Wen (2023).

### scripts

- `1.download.py`: script to download the waveform data
- `2.preprocess.py`: script to perform data preprocessing
- `3.pickphase.py`: pick phases manually
- `4.Fig1.py`: script for data analysis and visualization
- `fresnelzone.py`: simple script to calculate the Fresnel zone
- `helpers.py`: some helper functions used in other scripts

### Data files

- `waveforms`: waveform data (in SAC format) used in this study
- `Fig1.pdf`: Fig. 1 in this study

## Related References

- **Original Paper**:
  Yang, Y., & Song, X. (2022).
  Inner Core Rotation Captured by Earthquake Doublets and Twin Stations.
  *Geophysical Research Letters*, 49(12), e2022GL098393.
  https://doi.org/10.1029/2022GL098393
- **Comment**:
  Tian, D., & Wen, L. (2023).
  Comment on "Inner Core Rotation Captured by Earthquake Doublets and Twin Stations" by Yang and Song.
  *Geophysical Research Letters*, 50(15), e2023GL103173.
  https://doi.org/10.1029/2023GL103173
- **Reply**:
  Yang, Y., & Song, X. (2023).
  Misinterpreted seismic evidence for localized rapid changes of the inner core boundary surface.
  *Geophysical Research Letters*, 50(15), e2023GL104728.
  https://doi.org/10.1029/2023GL104728
