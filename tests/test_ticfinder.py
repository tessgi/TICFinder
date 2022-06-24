import numpy as np
import pandas as pd
import pytest

from ticfinder import TICFinder, __version__


def test_version():
    assert __version__ == "0.1.0"


@pytest.mark.remote_data
def test_ticfinder():
    tf = TICFinder(ra=84.2911880010838, dec=-80.4691198186792, magnitude=5.5)
    tf.get_tics()
    assert tf.tic[0] == 261136679
    df = pd.DataFrame(
        np.hstack([84.2911880010838, -80.4691198186792, 5.5])[None, :],
        columns=["ra", "dec", "tmag"],
    )
    tf = TICFinder.from_pandas(df)
    tf.get_tics()
    assert tf.tic[0] == 261136679
    tf.to_csv()
    tf = TICFinder(ra=84.2911880010838, dec=-80.4691198186792)
    tf.get_tics()
    assert tf.tic[0] == 261136679
    matches = tf.get_target_matches(0)
    assert len(matches) > 1
