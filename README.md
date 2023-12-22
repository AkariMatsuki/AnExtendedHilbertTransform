# AnExtendedHilbertTransform
Phase reconstruction from oscillatory data via an extended Hilbert transform (HT) method ([Matsuki, Kori and Kobayashi, Sci. Rep., 2023](https://www.nature.com/articles/s41598-023-30405-5)).  

Compared to the conventional HT method, the extended HT method can reconstruct phase more accurately from a sinusoidal signal $x(t) = A_L \cos (\hat{\omega} t + u(t))$, where $A_L$ is a slow amplitude and $u(t)$ is a small phase-modulation.


## Requirement
* Python 3.8.8
* NumPy 1.20.1
* SciPy 1.6.2
* Matplotlib 3.3.4
* sdeint 0.2.2

## Demo
See [Demo_An_Extended_Hilbert_Transform.ipynb](https://github.com/treepineakari1104/AnExtendedHilbertTransform/blob/main/Demo_An_Extended_Hilbert_Transform.ipynb)

## Usage
Reconstruct the phase via the conventional and the extended HT method by
```
python3 main.py "[file name]" tau
```

* file name: Name of the file(.txt) of the observed signal.
* tau: Sampling interval of the observed signal (float).

For example, 
```
python3 main.py "x_QuasiPeriodic.txt" 0.01
```
Two files will be generated:
* phase_exth.txt : phase reconstructed via the extended HT method.
* phase_h.txt : phase reconstructed via the conventional HT method.

## Author
Akari Matsuki
mail to : akarimatsuki114[at]gmail.com
