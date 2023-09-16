# PHAIN (PHase-aware Audio INpainter)
**Tomoro Tanaka (Department of Intermedia Art and Science, Waseda University, Tokyo, Japan)**\

This README file describes the MATLAB codes provided to test, analyze, and evaluate our proposed method, PHAINs.\
PHAIN is an audio declipping method introduced in the following paper
>[1] Tomoro Tanaka, Kohei Yatabe, and Yasuhiro Oikawa, "PHAIN: Audio Inpainting via Phase-aware Optimization with Instantaneous Frequency".

## Requirements
The codes were developed in MATLAB version R2023a and have been tested in R2023a and R2022b.\
Some functions rely on 

1. MathWorks Toolbox: You are kindly requested to download some of them, such as 'Signal Processing Toolbox'.

2. Toolbox available online: This is available online under the [MIT license](https://opensource.org/licenses/mit-license.php).

- Instantaneous phase correction of complex spectrogram\
  This folder contains a set of MATLAB codes written as a supplementary material of the following tutorial paper explaining phase-related property of spectrograms:
  I already installed it so you can easily execute the codes. Plaese refer to https://doi.org/10/c3qb.

  >[2] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa, "Representation of complex spectrogram via phase conversion," Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

## Data
There is a folder to contain test data. You can download and save any files there.\


## Usage
Execute `main.m` to perform PHAINs.

- `PHAIN`
  - `PHAINmain.m`.
  - `CP.m` is the Chambolle--Pock algorithm.

- `utils`
  - `proj_Gamma.m`.
  - `shortenForDGT.m`.

- `phase_correction`


## License
See the file named `LICENSE.pdf`.
