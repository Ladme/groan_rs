## Performance of `groan_rs` compiled with `molly` vs compiled with `xdrfile`

Performed using `groan_rs v0.10-dev.1`.

The benchmarking (atomistic) system was composed of a small peptide (~400 atoms), a phospholipid membrane (~40,000 atoms), and water with ions (~54,300 atoms). Total size: ~94,700 atoms.
The trajectory was 5000 frames long (0-500 ns).

The program iterated through a trajectory and printed out the position of the first atom for each analyzed frame.
Run time was measured with `time` (so only an estimate of run time is provided). Benchmarked on Debian 12 with an 8-core Intel Core i7-10700 CPU and SK Hynix BC511 HFM512GDJTNI SSD.

The reported run time includes the reading of a tpr and an ndx file.

|   atom selection   |       iteration type      | run time (groan + molly) [s] | run time (groan + xdrfile) [s] | rel. run time (molly/xdrfile) |
|:------------------:|:-------------------------:|:----------------------------:|:------------------------------:|:-----------------------------:|
|         all        |         all frames        |             10.1             |              16.3              |              62%              |
|         all        |      every 5th frame      |              2.1             |               3.2              |              66%              |
|         all        | frames between 300-400 ns |              2.1             |               3.2              |              66%              |

|   atom selection   |       iteration type      | run time (groan + molly) [s] | run time (groan + xdrfile) [s] | rel. run time (molly/xdrfile) |
|:------------------:|:-------------------------:|:----------------------------:|:------------------------------:|:-----------------------------:|
| peptide + membrane |         all frames        |              4.0             |               ---              |              ---              |
| peptide + membrane |      every 5th frame      |              0.9             |               ---              |              ---              |
| peptide + membrane | frames between 300-400 ns |              0.9             |               ---              |              ---              |

|   atom selection   |       iteration type      | run time (groan + molly) [s] | run time (groan + xdrfile) [s] | rel. run time (molly/xdrfile) |
|:------------------:|:-------------------------:|:----------------------------:|:------------------------------:|:-----------------------------:|
|       peptide      |         all frames        |              0.3             |               ---              |              ---              |
|       peptide      |      every 5th frame      |              0.1             |               ---              |              ---              |
|       peptide      | frames between 300-400 ns |              0.1             |               ---              |              ---              |

|   atom selection   |       iteration type      | run time (groan + molly) [s] | run time (groan + xdrfile) [s] | rel. run time (molly/xdrfile) |
|:------------------:|:-------------------------:|:----------------------------:|:------------------------------:|:-----------------------------:|
|        water       |         all frames        |              8.7             |               ---              |              ---              |
|        water       |      every 5th frame      |              1.8             |               ---              |              ---              |
|        water       | frames between 300-400 ns |              1.8             |               ---              |              ---              |

Note that when water is selected for reading, almost the entire XTC frames have to be read since water atoms are positioned near the end of the system. Yet, the run time is lower due to the lower amount of data copying performed inside the `groan` library.