# AnExtendedHilbertTransform
Phase reconstruction from oscillatory data via an extended Hilbert transform method.

## Requirement
* Python 3
* NumPy
* SciPy
* Matplotlib
* sdeint

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
mail to : treepineakari[at]g.ecc.u-tokyo.ac.jp
