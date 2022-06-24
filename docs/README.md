<a href="https://github.com/tessgi/ticfinder/actions/workflows/tests.yml"><img src="https://github.com/tessgi/ticfinder/workflows/pytest/badge.svg" alt="Test status"/></a> [![Generic badge](https://img.shields.io/badge/documentation-live-blue.svg)](https://tessgi.github.io/ticfinder)

# TICFinder

This tool helps you go from J2000 RA, Dec and [optionally] magnitude of sources and convert to TIC numbers, using the most up-to-date TIC catalog.

# Example Usage

Using the `TICFinder` class, you can pass RA, Dec and optionally magnitude to search for TIC IDs. (**Note: You must use J2000 epoch RA and Dec**.)

To use the class, pass in your RA, Dec and optional magnitude, and then use the `get_tics` class method to query for the TIC IDs. You can then access the `tic` attribute, or use the `to_pandas` method to create a csv file of your input with TIC IDs.

```python
from ticfinder import TICFinder
tf = TICFinder(ra=84.2911880010838, dec=-80.4691198186792, magnitude=5.5).get_tics()
print(tf.tic)
>>> [261136679]
```

You can also load a dataframe or a csv file.

```python
tf = TICFinder.from_pandas(df)
tf.get_tics()
print(tf.tic)
>>> [261136679, 261136641, 261136690, ...]
tf.to_pandas()
```

`tf.to_pandas` will return a `pandas.DataFrame` with columns

```python
['TIC', 'RA', 'Dec', 'Tmag', 'input_magnitude', 'pix_sep', 'motion_from_2000_to_YEAR_in_pixels']
```
